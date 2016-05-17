#! /usr/bin/env python
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import copy_reg
import os
import pysam
import types
import sys
import numpy as np
from smart_open import smart_open
from optparse import OptionParser
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from Rhea_Chip.lib import cnv, DescribeArray
from Rhea_Chip.lib.cnv import _chrom_valued, join_numbers, join_ranges, SaveLoad, CNVdata
from Rhea_Chip.lib.cnv.nbinom_fit import NegativeBinomial


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
	cpu_n = min(int(mem / use_mem), cpu_n, cpu_limit) if cpu_limit else min(int(mem / use_mem), cpu_n)
	return max(cpu_n, 1)


def _pickle_method(m):
	if m.im_self is None:
		return getattr, (m.im_class, m.im_func.func_name)
	else:
		return getattr, (m.im_self, m.im_func.func_name)


copy_reg.pickle(types.MethodType, _pickle_method)


class WinCorrect(object):
	def __init__(self, **kwargs):
		self.LowDepCut = float(kwargs["low_dep_cut"])
		self.CorrectWinLen = int(kwargs["correct_win_len"])
		self.CorrectShiftLen = int(kwargs["correct_shift_len"])
		chroms = SaveLoad(os.path.abspath(kwargs["chromstat"]))
		s_chrom = str(kwargs["chrom"]).split(",") if kwargs["chrom"] else None
		self.chrom_stat = chroms.load()
		self.samples = sorted(self.chrom_stat.keys())
		self.contigs = sorted(self.chrom_stat[self.samples[0]].keys(), key=lambda x: _chrom_valued(x))
		if s_chrom is not None and len(s_chrom):
			self.contigs = filter(lambda x: x in s_chrom, self.contigs)
		self.indir = os.path.abspath(kwargs["indir"])
		self.bed = dict()
		regions = defaultdict(list)
		bed = open(os.path.abspath(kwargs["region"]), 'r')
		for line in bed:
			if line.startswith("#"):
				continue
			rows = line.strip().split("\t")
			if len(rows) < 3:
				continue
			chrom = str(rows[0])
			if chrom not in self.contigs:
				continue
			start = int(rows[1])
			stop = int(rows[2])
			regions[chrom].extend(range(start, stop + 1))
		for chrom in self.contigs:
			self.bed[chrom] = sorted(regions[chrom])

	def run(self, debug=False):
		pool = Pool(processes=cpu(use_mem=3221225472, cpu_limit=len(self.contigs)))
		for chrom in self.contigs:
			pool.apply_async(self.win_correct, args=(chrom,))
		pool.close()
		pool.join()
		for chrom in self.contigs:
			cnvdata = dict()
			nbarg = os.path.join(self.indir, "%s.nbinom.arg" % chrom)
			if os.path.isfile(nbarg):
				f_in = smart_open(nbarg)
				trials, best_probability, devi = f_in.readline().strip().split("\t")
				f_in.close()
				if not debug:
					os.remove(nbarg)
			else:
				continue
			for sample in self.samples:
				cnvdata[sample] = CNVdata()
				cnvdata[sample].trials = int(trials)
				cnvdata[sample].best_probability = float(best_probability)
				cnvdata[sample].min_devi = float(devi)
				cnvdata[sample].ploid = int(self.chrom_stat[sample][chrom].ploid)
				cnvdata[sample].regions = list()
				cnvdata[sample].data = list()
				dep_data = os.path.join(self.indir, sample,
				                        "%s.W%dS%d.fixdep.gz" % (chrom, self.CorrectWinLen, self.CorrectShiftLen))
				if not os.path.isfile(dep_data):
					continue
				with smart_open(dep_data) as f_in:
					for line in f_in:
						if line.startswith("#"):
							continue
						chrom, start, stop, deps = line.strip().split("\t")
						start = int(start)
						stop = int(stop)
						deps = int(deps)
						cnvdata[sample].regions.append([chrom, start, stop])
						cnvdata[sample].data.append(deps)
			c_stat = SaveLoad(os.path.join(self.indir, "%s.cnv.args" % chrom))
			c_stat.save(cnvdata)

	def win_correct(self, chrom):
		pos_filter = defaultdict(int)
		filterSampleInChrom = set()
		bed_chrom = self.bed[chrom]
		depth_dict = dict()
		regions = list()
		sample_dep = defaultdict(list)
		win_sift_dep = defaultdict(list)
		wsdep = defaultdict(list)
		samples_filter = filter(lambda i: self.chrom_stat[i][chrom].ploid > 0, self.samples)
		if len(samples_filter) < 2:
			return
		nbinom_data = list()
		nbinom_out = open(os.path.join(self.indir, "%s.nbinom.arg" % chrom), 'w')
		chr_cor = open(os.path.join(self.indir, "%s_W%dS%d.cor" % (chrom, self.CorrectWinLen, self.CorrectShiftLen)),
		               'w')
		ws_deps = {sample: os.path.join(self.indir, sample, "%s.W%dS%d.fixdep" %
		                                (chrom, self.CorrectWinLen, self.CorrectShiftLen))
		           for sample in samples_filter}
		ws_dep = {sample: open(ws_deps[sample], 'w') for sample in samples_filter}
		for sample in samples_filter:
			chrom_d = list()
			fixdeps = pysam.TabixFile(os.path.join(self.indir, "{0}/{0}.Fixdep.tsv.gz".format(sample)))
			depth_dict[sample] = fixdeps
			for line in fixdeps.fetch(chrom):
				rows = line.strip().split("\t")
				pos = int(rows[1])
				c_d = int(rows[-1])
				chrom_d.append([pos, c_d])
			depths = DescribeArray(chrom_d, col=1)
			for (pos, dep) in chrom_d:
				if dep < 0.6 * self.LowDepCut * depths.average:
					pos_filter[pos] += 2
				elif dep < self.LowDepCut * depths.average:
					pos_filter[pos] += 1
		for pos, number in pos_filter.iteritems():
			if number >= len(samples_filter):
				try:
					bed_chrom.remove(pos)
				except ValueError:
					continue
		bed_chrom = list(join_ranges(join_numbers(bed_chrom), offset=self.CorrectWinLen))
		for s, e in bed_chrom:
			if e - s < self.CorrectShiftLen:
				continue
			for win in xrange(s, e, self.CorrectShiftLen):
				end_p = min(win + self.CorrectWinLen - 1, e)
				regions.append((chrom, win, end_p))
				for sample in samples_filter:
					depth = depth_dict[sample]
					lines = list(depth.fetch(chrom, win, end_p)) or ["-1\t-1\t0\n"]
					mdep = sum([int(line.strip().split("\t")[-1]) for line in lines]) / float(len(lines))
					mdep = round(mdep, 2)
					sample_dep[sample].append(mdep)
					win_sift_dep[(chrom, win, end_p)].append(mdep)
		mdep_chr = {i: sum(j) / len(j) for i, j in sample_dep.iteritems() if len(j)}
		cormtx = np.corrcoef([sample_dep[i] for i in samples_filter])
		chr_cor.write(chrom + "\t" + "\t".join(samples_filter) + '\n')
		for sn in range(len(samples_filter)):
			chr_cor.write("\t".join([samples_filter[sn]] + map(str, cormtx[sn])))
			mcor = (cormtx[sn].sum() - 1.0) / (len(samples_filter) - 1.0)
			if mcor < 0.6:
				filterSampleInChrom.add(samples_filter[sn])
				chr_cor.write("\tLow correlation\n")
			else:
				chr_cor.write("\n")
		chr_cor.write("\n")
		if len(filterSampleInChrom) / float(len(samples_filter)) > 0.6:
			return
		for r in regions:
			chrom, start, stop = r
			d = win_sift_dep[r]
			dr = list()
			ndr = list()
			for sn in range(len(samples_filter)):
				sample = samples_filter[sn]
				if sample in filterSampleInChrom:
					continue
				contorl = self.LowDepCut * mdep_chr[sample]
				dr.append(d[sn] / mdep_chr[sample]) if mdep_chr[sample] > 0 else 0
				if d[sn] > contorl:
					ndr.append(d[sn] / mdep_chr[sample])
			mdr = np.median(np.array(ndr)) if len(ndr) > 3 else 1.0
			if mdr < self.LowDepCut:
				mdr = 1.0
			for sn in range(len(samples_filter)):
				d[sn] /= mdr
				wsdep[samples_filter[sn]].append(d[sn])
		for s, d in wsdep.iteritems():
			if s in filterSampleInChrom:
				continue
			tmp_des = DescribeArray(d)
			d = [min(5.0 * tmp_des.average, i) for i in d]
			tmp_des = DescribeArray(d)
			d = [int(i * 60.0 / tmp_des.average + 0.5) for i in d]
			if len(regions) == len(d):
				for i in range(len(regions)):
					reg = "\t".join(map(str, regions[i]))
					ws_dep[s].write("%s\t%i\n" % (reg, d[i]))
					nbinom_data.append(d[i])
		chr_cor.close()
		for s, d in ws_deps.iteritems():
			ws_dep[s].close()
			depth_dict[s].close()
			if os.path.isfile(d):
				_ = pysam.tabix_index(d, seq_col=0, start_col=1, end_col=2, force=True)
		nbarg = NegativeBinomial(nbinom_data)
		nbarg.nbinom_fit()
		trials = nbarg.trials
		best_probability = nbarg.best_probability
		min_devi = nbarg.mindevi
		nbinom_out.write("\t".join(map(str, [trials, best_probability, min_devi])) + '\n')
		nbinom_out.close()


def main():
	usage = 'Usage: %prog [-h] [--version] -i [fixdep dir] -c [chrom stat file] [options]'
	parser = OptionParser(usage=usage, version=cnv.__version__, description=cnv.__description__, epilog=cnv.__author__)
	parser.add_option('--low_dep_cut', dest='low_dep_cut', help='low depth cut off, default: 0.2', default=0.2)
	parser.add_option('--correct_win_len', dest='correct_win_len', help='windows length, default: 30', default=30)
	parser.add_option('--correct_shift_len', dest='correct_shift_len', help='Shift lenght, default: 25', default=25)
	parser.add_option('-c', '--chromstat', metavar='FILE', dest='chromstat', help='chrom stat file')
	parser.add_option('-i', '--indir', metavar='DIR', dest='indir', help='The fixdep dir')
	parser.add_option('-b', '--bed', metavar='FILE', dest='region', help='The bed region')
	parser.add_option('-C', '--chrom', dest='chrom', help='special chromosome, default: all', default=None)
	(options, args) = parser.parse_args()
	if not options.indir or not options.chromstat or not options.region:
		parser.print_help()
		return 'Expected arguments lost !!!'
	jobs = WinCorrect(**options.__dict__)
	jobs.run()


if __name__ == '__main__':
	sys.exit(main())
