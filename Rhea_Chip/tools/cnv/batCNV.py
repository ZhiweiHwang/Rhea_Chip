#! /usr/bin/env python
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
import pysam
import sys
import numpy as np
from optparse import OptionParser
from glob import glob
from multiprocessing import Pool, cpu_count
from smart_open import smart_open
from itertools import groupby
from collections import defaultdict
from Rhea_Chip.lib.cnv import _chrom_valued, __version__, __description__, __author__, SaveLoad, count_gc, join_ranges
from Rhea_Chip.lib.annovar import AnnoVar
from Rhea_Chip.lib.cnv.nbinom_fit import NegativeBinomial
from Rhea_Chip.lib.cnv.hmm import ViterbiTraining
from Rhea_Chip.lib.hgvs import HGVS
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

_AnnotationDB = defaultdict(set)


def cpu(use_mem=1073741824, cpu_limit=0):
	mem = sys.maxint
	if os.path.exists("/proc/meminfo"):
		f = open("/proc/meminfo")
		for line in f:
			if len(line) < 2:
				continue
			name, var = line.split(':')[0:2]
			if name == 'MemFree':
				mem = float(var.split()[0]) * 1024.0
				break
		f.close()
	cpu_n = cpu_count() - 2
	cpu_n = min(int(mem / use_mem - 0.5), cpu_n, cpu_limit) if cpu_limit else min(int(mem / use_mem - 0.5), cpu_n)
	return max(cpu_n, 1)


def hmm_cnv(dep_data, regions, best_probability, trials, ploid=2.0, output=None, contral_wins=5):
	output = os.path.abspath(output) if output is not None else sys.stdout
	final_cnv = smart_open(output, 'w')
	cne = ViterbiTraining(dep_data, ploid, best_probability, trials)
	copy_est = cne.train()
	if np.count_nonzero(copy_est != ploid) / float(len(copy_est)) < 0.1 and len(copy_est) > 800:
		iterations = 0
		lastDif = differences = len(copy_est)
		n_copy_est = np.copy(copy_est)
		tmp_arg = [[0, 0, 0], [0, 0, 0]]
		while differences > 0 and iterations < 100:
			ndep = [int(ploid * dep_data[i] / float(copy_est[i]) + 0.5) for i in range(len(copy_est)) if copy_est[i]]
			nbarg = NegativeBinomial(ndep)
			nbarg.nbinom_fit()
			tmp_b = float(nbarg.best_probability)
			tmp_t = int(nbarg.trials)
			devi = float(nbarg.mindevi)
			cne = ViterbiTraining(dep_data, ploid, tmp_b, tmp_t)
			copy_est = cne.train(copy_est)
			iterations += 1
			differences = np.count_nonzero(np.array(n_copy_est) != copy_est)
			n_copy_est = np.copy(copy_est)
			if differences == lastDif:
				if (tmp_arg[0][2] == devi) and (tmp_arg[0][1] == tmp_t) and (tmp_arg[0][0] == tmp_b):
					break
			lastDif = differences
			tmp_arg = [[tmp_arg[1][0], tmp_arg[1][1], tmp_arg[1][2]], [tmp_b, tmp_t, devi]]
	trans_p = cne.mlEstimate(copy_est)
	cnvs = cne.posterior_decoding(trans_p, ploid)
	chrom = str(regions[0][0])
	for k, g in groupby(zip(regions, cnvs), lambda x: x[1][1]):
		if k == ploid:
			continue
		elems = list(g)
		for start, stop in join_ranges([d[0][1:] for d in elems]):
			p = [float(d[1][0]) for d in elems if d[0][1] >= start and d[0][2] <= stop]
			mpp = round(sum(p) / len(p), 3)
			if mpp < 0.95 or len(p) < contral_wins:
				continue
			bp_len = stop - start + 1
			mut_type = "gain" if k > ploid else "loss"
			final_cnv.write("\t".join(map(str, [chrom, ploid, start, stop, bp_len, k, mut_type, mpp])) + '\n')
	final_cnv.close()


class CNVAnnotation(object):
	def __init__(self, reference, annotationdb=_AnnotationDB):
		self.refer = pysam.FastaFile(os.path.abspath(reference))
		self._dbhandles = list()
		dbtitle = defaultdict(float)
		for dbconfigset in annotationdb.values():
			for dbconfig in dbconfigset:
				dbs = AnnoVar(dbconfig)
				if "Start" in dbs.search_value and "Stop" in dbs.search_value:
					v = filter(lambda x: x not in ["Start", "Stop"], dbs.search_value.values())
					v.append("_".join([str(dbs.db), "RegionCover"]))
					v.sort()
				else:
					v = sorted(dbs.search_value.values())
				t = "\t".join(map(str, v))
				dbtitle[t] = dbs.order
				self._dbhandles.append(dbs)
		self.dbtitle = "\t".join([j for j in sorted(dbtitle.keys(), key=lambda x: dbtitle[x])])

	def __del__(self):
		self.refer.close()
		for dbs in self._dbhandles:
			dbs.__del__()

	def close(self):
		self.__del__()

	def dbanno(self, variation):
		dbinfo = defaultdict(set)
		for dbs in self._dbhandles:
			try:
				infos = dbs.fetch(**variation)
				for info in infos:
					if "Start" in dbs.search_value and "Stop" in dbs.search_value:
						r_s = int(info.Start)
						r_e = int(info.Stop)
						radio = (min(r_e, variation["Stop"]) - max(variation["Start"], r_s) + 1) / float(r_e - r_s + 1)
						if radio < 0.001:
							continue
						dbinfo["_".join([dbs.db, "RegionCover"])].add(
							"%s:%i-%i[%.4f]" % (variation["Chrom"], r_s, r_e, radio))
					for k, v in info.__dict__.iteritems():
						if k == "Start" or k == "Stop":
							continue
						if v != ".":
							dbinfo[k].add(v)
			except Exception:
				continue
		return dbinfo


class CNVAnalysis(object):
	global _AnnotationDB

	def __init__(self, **kwargs):
		self.outdir = os.path.abspath(kwargs["indir"])
		self.win_len = int(kwargs["correct_win_len"]) or 30
		self.shift_len = int(kwargs["correct_shift_len"]) or 25
		self.contral_wins = int(100.0 / self.shift_len + 0.5) + 1
		chroms = str(kwargs["chrom"]).split(",") if kwargs["chrom"] else None
		samples = str(kwargs["sample"]).split(",") if kwargs["sample"] else None
		all_samples = set()
		contigs = list()
		self.cnvdata = defaultdict(dict)
		self.sample_win_data = defaultdict(dict)
		for cnv_data in glob(os.path.join(self.outdir, "chr*.cnv.args")):
			chrom = ".".join(os.path.basename(cnv_data).split(".")[0:-2])
			if chroms is not None and chrom not in chroms:
				continue
			cnvdata = SaveLoad(cnv_data)
			cnvdata = cnvdata.load()
			contigs.append(chrom)
			for sample in cnvdata.keys():
				dep_f = os.path.join(self.outdir, sample, "%s.W%iS%i.fixdep.gz" % (chrom, self.win_len, self.shift_len))
				if os.path.isfile(dep_f) and os.path.isfile(dep_f + '.tbi'):
					self.sample_win_data[chrom][sample] = dep_f
				if samples is not None and sample not in samples:
					continue
				all_samples.add(sample)
				self.cnvdata[sample][chrom] = cnvdata[sample]
		self.samples = sorted(all_samples)
		self.contigs = sorted(contigs, key=lambda x: _chrom_valued(x))
		databases = os.path.abspath(kwargs["dbdir"])
		t_db = os.path.abspath(kwargs["transdb"]) if "transdb" in kwargs else os.path.join(databases, "transdb",
		                                                                                   "ncbi_anno_rel104.dbref.db")
		for db in glob(os.path.join(databases, "*", "*.cnvdb.config")):
			db = os.path.abspath(db)
			dbname = os.path.basename(os.path.dirname(db))
			_AnnotationDB[dbname].add(db)
		self.reference = os.path.abspath(kwargs["reference"]) if kwargs["reference"] else \
			os.path.join(databases, 'aln_db/hg19/hg19_chM_male_mask.fa')
		self.DBAnno = CNVAnnotation(self.reference, _AnnotationDB)
		self.HGVS = HGVS(t_db)

	def call_cnvs(self, sample):
		pool = Pool(processes=cpu(use_mem=2147483648, cpu_limit=len(self.contigs)))
		for chrom in self.contigs:
			dep_data = self.cnvdata[sample][chrom].data
			ploid = self.cnvdata[sample][chrom].ploid
			if ploid < 1 or not len(dep_data):
				continue
			best_p = self.cnvdata[sample][chrom].best_probability
			trials = self.cnvdata[sample][chrom].trials
			regions = self.cnvdata[sample][chrom].regions
			output = os.path.join(self.outdir, sample, "%s.cnv" % chrom)
			hmm_cnv(dep_data, regions, best_p, trials, ploid, output, self.contral_wins)
			pool.apply_async(hmm_cnv, (dep_data, regions, best_p, trials, ploid, output, self.contral_wins))
		pool.close()
		pool.join()

	def z_score(self, chrom, start, stop, sample):
		dep_data = list()
		deps = list()
		dep_handle = pysam.Tabixfile(self.sample_win_data[chrom][sample])
		for i in dep_handle.fetch(chrom):
			p_s, p_e, p_d = map(int, i.strip().split("\t")[1:4])
			dep_data.append(p_d)
			if (p_s <= start <= p_e) or (start <= p_s <= stop) or (p_s <= stop <= p_e):
				deps.append(p_d)
		dep_handle.close()
		if len(deps) == 0:
			return 0
		d_m = np.mean(dep_data)
		s_d = np.std(dep_data)
		score = 0.0
		for i in deps:
			score += (i - d_m + 0.0) / s_d
		return round(score / len(deps), 4)

	def plot(self, sample, f_in):
		fin = smart_open(f_in)
		f_out = os.path.join(os.path.dirname(f_in), sample + ".cnv.out.pdf")
		pdf = PdfPages(f_out)
		for line in fin.readlines():
			rows = line.strip().split("\t")
			try:
				chrom = str(rows[0])
				start = max(int(rows[2]) - self.win_len, 0)
				stop = int(rows[3]) + self.win_len
			except ValueError:
				continue
			plt.figure()
			samples = sorted(self.sample_win_data[chrom].keys())
			deps = list()
			for s in samples:
				d = list()
				dep_handle = pysam.Tabixfile(self.sample_win_data[chrom][s])
				for lines in dep_handle.fetch(chrom, start, stop):
					p_s, p_e, p_d = map(int, lines.strip().split("\t")[1:4])
					d.append(p_d)
				dep_handle.close()
				deps.append(d)
			deps = np.array(deps)
			deps = np.log2(deps / deps.mean(axis=0))
			for s in range(len(samples)):
				d = deps[s]
				if samples[s] == sample:
					plt.plot(d, color="k")
				else:
					plt.plot(d, color="m")
			plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
			plt.title("%s_%s_%s_%s" % (rows[0], rows[2], rows[3], rows[5]))
			pdf.savefig()
			plt.close()
		pdf.close()

	def annotation(self, sample, debug=False):
		output = os.path.join(self.outdir, sample, "%s.cnv.anno.tsv" % sample)
		f_out = smart_open(output, 'w')
		titles = ["#Chrom", "Ploid", "Start", "Stop", "length", "copyNumber", "Mtype", "meanP", "Z-score", "GCpr",
		          "MutationName"]
		dbtitle = self.DBAnno.dbtitle.split("\t")
		titles.extend(dbtitle)
		f_out.write("\t".join(titles) + '\n')
		for chrom in self.contigs:
			cnvs = os.path.join(self.outdir, sample, "%s.cnv" % chrom)
			if not os.path.exists(cnvs):
				continue
			f_in = smart_open(cnvs)
			for line in f_in:
				rows = line.strip().split("\t")
				try:
					chrom = str(rows[0])
					start = int(rows[2])
					stop = int(rows[3])
					mtype = str(rows[6])
				except ValueError:
					continue
				z_s = self.z_score(chrom, start, stop, sample)
				gcr = count_gc(self.DBAnno.refer.fetch(chrom, start, stop))
				if not (0.3 <= gcr <= 0.7):
					continue
				variation = dict([("Chrom", chrom), ("Start", start), ("Stop", stop), ("Mtype", mtype)])
				m_name = set()
				for hgvs in self.HGVS.annobed(chrom, start, stop):
					trans = str(hgvs.Transcript)
					gene = str(hgvs.geneSym)
					chgvs = str(hgvs.cHgvs)
					protein = str(hgvs.Protein)
					exons = str(hgvs.ExonRegions)
					mess = ":".join(filter(lambda x: x != ".", [trans, protein, gene, chgvs, exons]))
					m_name.add(mess)
				dbinfo = self.DBAnno.dbanno(variation)
				anno_message = [str(i) for i in rows]
				anno_message.append(str(z_s))
				anno_message.append(str(gcr))
				anno_message.append("|".join(m_name))
				for i in dbtitle:
					if i in dbinfo:
						anno_message.append("|".join(dbinfo[i]))
					else:
						anno_message.append(".")
				f_out.write("\t".join(anno_message) + '\n')
			f_in.close()
			if not debug:
				os.remove(cnvs)
		f_out.close()
		return output


def main():
	usage = 'Usage: %prog [-h] [--version] -i [fixdep dir] -c [chrom stat file] [options]'
	parser = OptionParser(usage=usage, version=__version__, description=__description__, epilog=__author__)
	parser.add_option('-r', '--reference', dest='reference', help='human reference', default=None)
	parser.add_option('-i', '--indir', metavar='DIR', dest='indir', help='The fixdep dir')
	parser.add_option('-d', '--db', metavar='DIR', dest='dbdir', help='The annotation database dirs')
	parser.add_option('-c', '--chrom', dest='chrom', help='special chromosome, default: all', default=None)
	parser.add_option('-s', '--sample', dest='sample', help='special sample, default: all', default=None)
	parser.add_option('--correct_win_len', dest='correct_win_len', help='windows length, default: 30', default=30)
	parser.add_option('--correct_shift_len', dest='correct_shift_len', help='Shift lenght, default: 25', default=25)
	parser.add_option('-p', '--plot', help='Plot exon depth graph', dest='plot', default=False, action="store_true")
	(options, args) = parser.parse_args()
	if not options.indir or not options.dbdir:
		parser.print_help()
		return 'Expected arguments lost !!!'
	jobs = CNVAnalysis(**options.__dict__)
	plot = options.plot
	for sample in jobs.samples:
		jobs.call_cnvs(sample)
		output = jobs.annotation(sample)
		if plot:
			jobs.plot(sample, output)

if __name__ == '__main__':
	sys.exit(main())
