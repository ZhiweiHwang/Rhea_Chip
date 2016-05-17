#! /usr/bin/env python
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import pysam
import os
import sys
import numpy as np
import statsmodels.api as sm
from smart_open import smart_open
from Rhea_Chip.lib import cnv
from Rhea_Chip.lib import DescribeArray
from Rhea_Chip.lib.cnv import _chrom_valued, SaveLoad
from optparse import OptionParser

lowess = sm.nonparametric.lowess


def unique_rows(a):
	a = np.ascontiguousarray(a)
	unique_a = np.unique(a.view([('', a.dtype)] * a.shape[1]))
	return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))


def gc_correct(**kwargs):
	depthf = os.path.abspath(kwargs["input"])
	if not os.path.isfile(depthf) or not os.path.isfile(depthf + '.tbi'):
		return
	sample = str(kwargs['sample']) if kwargs['sample'] else os.path.basename(depthf).split(".")[0]
	outdir = os.path.join(os.path.abspath(kwargs['outdir']), sample)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	out = os.path.join(outdir, "%s.Fixdep.tsv" % sample)
	wins = SaveLoad(os.path.abspath(kwargs["wingc"]))
	wingc = wins.load()
	poss = SaveLoad(os.path.abspath(kwargs["posgc"]))
	posgc = poss.load()
	chroms = SaveLoad(os.path.abspath(kwargs["chromstat"]))
	chrom_stat = chroms.load()
	f_out = smart_open(out, 'w')
	f_out.writelines("#Chrom\tPos\tFixDepth\n")
	gc_depth = list()
	dep_f = pysam.TabixFile(depthf)
	for rows, gc_content in sorted(wingc.iteritems(), key=lambda x: (_chrom_valued(x[0][0]), x[0][1])):
		chrom = rows[0]
		start = rows[1]
		stop = rows[2] - 1
		if chrom not in chrom_stat[sample] or chrom_stat[sample][chrom] < 1:
			continue
		try:
			depths = [int(line.strip().split("\t")[-1]) for line in dep_f.fetch(chrom, start, stop)]
			win_mean_dep = min(sum(depths) / float(len(depths)), 6.0 * chrom_stat[sample][chrom].average)
			win_mean_dep *= 2.0 / chrom_stat[sample][chrom].ploid
		except Exception:
			win_mean_dep = 0.0
		gc_depth.append([gc_content, win_mean_dep])
	gc_depth = DescribeArray(gc_depth, col=1)
	gcdep = gc_depth.array[gc_depth.array[:, 1] > 0.05 * gc_depth.median]
	prd = unique_rows(lowess(gcdep[:, 1], gcdep[:, 0], frac=0.25))
	mdp = np.median(prd[:, 1])
	if mdp <= 0.0:
		raise ValueError("Sample %s depth file Error !" % depthf)
	lgc = gcl = max(10000, int(prd[:, 0].max() * 10000))
	loe = [-0.0001, ] * gcl
	gcj = 0
	for gc, dp in prd:
		gcj = int(round(gc, 4) * 10000)
		if gcj < gcl:
			gcl = gcj
		loe[gcj] = mdp / float(dp) if dp > 0 else 1.0
	for gc in xrange(gcl):
		loe[gc] = min(loe[gcl], 10.0)
	for i in xrange(gcl + 1, gcj):
		if loe[i] < 0:
			ls = i - 1
			lv = loe[i - 1]
			rs = i + 1
			while loe[rs] < 0 and rs < len(loe):
				rs += 1
			rv = loe[rs]
			loe[i] = min((lv + (rs - float(ls)) * rv) / (rs - float(ls) + 1.0), 10.0)
	for i in xrange(gcj + 1, lgc):
		loe[i] = min(loe[gcj], 10.0)
	for line in dep_f.fetch():
		rows = line.strip().split("\t")
		chrom = str(rows[0])
		pos = int(rows[1])
		deps = int(rows[-1])
		try:
			fixdeps = int(deps * loe[int(round(posgc[chrom][pos], 4) * 10000)])
		except KeyError:
			continue
		f_out.writelines("\t".join(map(str, [chrom, pos, fixdeps])) + '\n')
	f_out.close()
	dep_f.close()
	_ = pysam.tabix_index(out, seq_col=0, start_col=1, end_col=1, force=True)


def main():
	usage = 'Usage: %prog [-h] [--version] --w [winGC file] -i [depth file] -p [posGC file] -s [sample name] [options]'
	parser = OptionParser(usage=usage, version=cnv.__version__, description=cnv.__description__, epilog=cnv.__author__)
	parser.add_option('-w', '--wingc', metavar='FILE', dest='wingc', help='windows GC file')
	parser.add_option('-p', '--posgc', metavar='FILE', dest='posgc', help='pos GC file')
	parser.add_option('-i', '--input', metavar='FILE', dest='input', help='input depth')
	parser.add_option('-c', '--chromstat', metavar='FILE', dest='chromstat', help='chrom stat file')
	parser.add_option('-s', '--samplename', dest='sample', help='sample name', default=None)
	parser.add_option('-o', '--outdir', metavar='DIR', dest='outdir', default=os.getcwd(),
	                  help='The output file [ default: stdout ]')
	(options, args) = parser.parse_args()
	if not options.wingc or not options.posgc or not options.input or not options.chromstat:
		parser.print_help()
		return 'Expected arguments lost !!!'
	gc_correct(**options.__dict__)


if __name__ == '__main__':
	sys.exit(main())
