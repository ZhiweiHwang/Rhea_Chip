#! /usr/bin/env python
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！
"""
Remove repeat regions, and count other regions' GC content
"""
import os
import pysam
import sys
import numpy as np
from collections import defaultdict
from optparse import OptionParser
from smart_open import smart_open
from Rhea_Chip.lib import cnv, DescribeArray
from Rhea_Chip.lib.cnv import join_ranges, SaveLoad, count_gc

pos_gc = defaultdict(float)
win_gc = dict()


class RegionAnalysis(object):
	def __init__(self, reference, rmtkdb):
		self.reference = pysam.FastaFile(reference)
		self.rmtk = pysam.TabixFile(rmtkdb)
		self.contigs = set()

	def __del__(self):
		self.reference.close()
		self.rmtk.close()

	def analysis(self, chrom, start, stop, win_len=200, sift_len=20):
		global pos_gc, win_gc
		chrom = str(chrom)
		tmp_start = start = int(start)
		stop = int(stop)
		flank_stop = min(int(stop) + win_len + 1, self.reference.get_reference_length(chrom))
		try:
			rmtk = list(self.rmtk.fetch(chrom, start, stop))
			bases = self.reference.fetch(chrom, start, flank_stop)
		except Exception as err:
			raise ValueError(err)
		self.contigs.add(chrom)
		for pos in xrange(start, stop + 1):
			flank_b = max(pos - win_len / 2, 0)
			flank_e = min(flank_b + win_len + 1, self.reference.get_reference_length(chrom))
			base_gc = count_gc(self.reference.fetch(chrom, flank_b, flank_e).upper())
			pos_gc[chrom][pos] = base_gc
		feback = list()
		if len(rmtk):
			for lines in rmtk:
				rows = lines.strip().split("\t")
				begin = int(rows[1])
				end = int(rows[2])
				if begin > start:
					feback.append([start, begin - 1])
				else:
					begin = start
				end = min(stop, end)
				if begin > end:
					continue
				start = end + 1
		else:
			feback.append([start, stop])
		for s, e in list(join_ranges(feback, offset=60)):
			if e - s < sift_len:
				continue
			for win in xrange(s, e, sift_len):
				offset = win - tmp_start
				seq = bases[offset:offset + win_len].upper()
				gc_radio = count_gc(seq)
				if 0.1 < gc_radio < 0.9:
					win_gc[(chrom, win, win + win_len - 1)] = gc_radio
				else:
					continue

	def chrom_stat(self, depth_files):
		smd = dict()
		gender = dict()
		ChromStat = defaultdict(dict)
		for dep in depth_files:
			if not os.path.isfile(dep) or not os.path.isfile(dep + '.tbi'):
				continue
			sample = os.path.basename(dep).split(".")[0]
			depth = pysam.TabixFile(dep)
			all_chrom_ave = list()
			for chrom in depth.contigs:
				if chrom not in self.contigs:
					continue
				depths = [float(line.strip().split("\t")[-1]) for line in depth.fetch(reference=chrom)]
				ChromStat[sample][chrom] = DescribeArray(depths)
				all_chrom_ave.append(ChromStat[sample][chrom].average)
			smd[sample] = np.median(np.array(all_chrom_ave)) if len(all_chrom_ave) else 0.0
			depth.close()
		for chrom in self.contigs:
			dr = list()
			for sample, sample_mean_depth in smd.iteritems():
				if chrom not in ChromStat[sample]:
					ChromStat[sample][chrom] = DescribeArray([0.0])
				c_s_c = float(ChromStat[sample][chrom].average)
				dr.append(c_s_c / smd[sample]) if smd[sample] else 0
			mdr = np.median(np.array(dr, dtype=object))
			for sample in smd.keys():
				c_s_c = float(ChromStat[sample][chrom].average)
				fxmed = c_s_c / mdr if mdr > 0.0 else 0.0
				if fxmed < 0.25 * smd[sample]:
					ChromStat[sample][chrom].ploid = 0.0
				elif fxmed < 0.55 * smd[sample]:
					ChromStat[sample][chrom].ploid = 1.0
				elif fxmed < 1.45 * smd[sample]:
					ChromStat[sample][chrom].ploid = 2.0
				elif fxmed < 1.95 * smd[sample]:
					ChromStat[sample][chrom].ploid = 3.0
				else:
					ChromStat[sample][chrom].ploid = min(round(2 * fxmed / smd[sample]), 5.0)
		if 'chrX' in self.contigs or 'chrY' in self.contigs:
			for dep in depth_files:
				if not os.path.isfile(dep) or not os.path.isfile(dep + '.tbi'):
					continue
				sample = os.path.basename(dep).split(".")[0]
				depth = pysam.TabixFile(dep)
				n_dep = np.median(np.array([ChromStat[sample][i].average for i in depth.contigs
				                            if i not in ["chrX", "chrY"]]))
				if 'chrX' in depth.contigs:
					x_dep = ChromStat[sample]["chrX"].average
					if 'chrY' not in depth.contigs:
						if abs(x_dep - n_dep) > abs(x_dep - 0.5 * n_dep):
							gender[sample] = "male"
							ChromStat[sample]["chrX"].ploid = 1.0
						else:
							gender[sample] = "female"
					else:
						y_dep = ChromStat[sample]["chrY"].average
						if y_dep > 0.15 * n_dep:
							gender[sample] = "male"
							ChromStat[sample]["chrX"].ploid = 2.0 if x_dep > 0.75 * n_dep else 1.0
							ChromStat[sample]["chrY"].ploid = 1.0
						else:
							gender[sample] = "female"
							ChromStat[sample]["chrX"].ploid = 1.0 if x_dep < 0.7 * n_dep else 2.0
							ChromStat[sample]["chrY"].ploid = 0.0
				else:
					y_dep = ChromStat[sample]["chrY"].average
					if y_dep > 0.1 * n_dep:
						gender[sample] = "male"
						ChromStat[sample]["chrY"].ploid = 1.0
					else:
						gender[sample] = "female"
						ChromStat[sample]["chrY"].ploid = 0.0
				depth.close()
		return ChromStat


def bedAnalysis(**kwargs):
	global pos_gc, win_gc
	bed = os.path.abspath(kwargs["bed"])
	reference = os.path.abspath(kwargs["reference"])
	db = os.path.abspath(kwargs["db"])
	outdir = os.path.abspath(kwargs["outdir"])
	winlen = int(kwargs["winlen"]) if "winlen" in kwargs else 200
	siftlen = int(kwargs["siftlen"]) if "siftlen" in kwargs else 20
	depth_f = [os.path.abspath(i) for i in kwargs["depthfile"].split(",") if os.path.isfile(i)]
	model = RegionAnalysis(reference, db)
	bed_gc_out = SaveLoad(os.path.join(outdir, "win.gc"))
	pos_gc_out = SaveLoad(os.path.join(outdir, "pos.gc"))
	chrom_stat = SaveLoad(os.path.join(outdir, "chrom.stat"))
	with smart_open(bed) as f_in:
		for line in f_in:
			rows = line.strip().split("\t")
			chrom = str(rows[0])
			if chrom not in pos_gc:
				pos_gc[chrom] = dict()
			start = int(rows[1])
			stop = int(rows[2]) + 1
			try:
				model.analysis(chrom, start, stop, winlen, siftlen)
			except ValueError:
				continue
	bed_gc_out.save(win_gc)
	pos_gc_out.save(pos_gc)
	chrom_stat.save(model.chrom_stat(depth_f))
	model.__del__()
	return bed_gc_out.fname, pos_gc_out.fname, chrom_stat.fname


def main():
	usage = 'Usage: %prog [-h] [--version] --bed [target bed] -r [ref seq] -d [RepeatMasker databases] [options]'
	parser = OptionParser(usage=usage, version=cnv.__version__, description=cnv.__description__, epilog=cnv.__author__)
	parser.add_option('-b', '--bed', metavar='FILE', dest='bed', help='bed file')
	parser.add_option('-r', '--reference', metavar='FILE', dest='reference', help='humman reference')
	parser.add_option('-d', '--db', dest='db', metavar='FILE', help='RepeatMasker databases')
	parser.add_option('-f', '--files', dest='depthfile', help='depth files, join by ","')
	parser.add_option('-o', '--out', dest='outdir', default=os.getcwd(),
	                  help='The output dir [ default: %s ]' % os.getcwd())
	parser.add_option('-w', '--winlen', dest='winlen', default=200, help='The windowns length [ default: 200 ]')
	parser.add_option('-s', '--siftlen', dest='siftlen', default=20, help='The windowns sift length [ default: 20 ]')
	(options, args) = parser.parse_args()
	if not options.bed or not options.reference or not options.db or not options.depthfile:
		parser.print_help()
		return 'Expected arguments lost !!!'
	bedAnalysis(**options.__dict__)


if __name__ == '__main__':
	sys.exit(main())
