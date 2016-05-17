#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# The MIT License (MIT)
# Copyright (c) 2015 Joey Hwong
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
# associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions: The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import os
import sys
from smart_open import smart_open
from optparse import OptionParser, OptionGroup
from Rhea_Chip.lib import cnv, BASH
from Rhea_Chip.lib.cnv.bedAnalysis import bedAnalysis
from Rhea_Chip.lib.cnv.GCcorrect import GCcorrect
from Rhea_Chip.lib.cnv.WinCorrect import WinCorrect
from Rhea_Chip.lib.cnv.batCNV import batCNV
from Rhea_Chip.lib.CreateFolder import FolderMaker

Rhea_Chip_Home = os.getenv('Rhea_Chip_Home') or os.path.dirname(os.path.abspath(__file__))


class Pipeline(object):
	"""
	create all scripts for analysis
	"""

	def __init__(self, **kwargs):
		self.region = os.path.abspath(kwargs["region"])
		sample_f = os.path.abspath(kwargs["samplelist"])
		depth_dir = os.path.abspath(kwargs["depths"]) if kwargs["depths"] is not None else os.path.join(
			os.path.dirname(sample_f), 'QC')
		self.reference = os.path.abspath(kwargs["reference"])
		self.rmsk = os.path.abspath(kwargs["repeatdb"])
		self.outdir = os.path.abspath(kwargs["outdir"])
		self.script = os.path.join(self.outdir, 'script')
		self.log = os.path.join(self.script, 'log')
		self.winlen = int(kwargs["winlen"]) or 30
		self.siftlen = int(kwargs["siftlen"]) or 25
		self.SamplesDepth = dict()
		self.contigs = set()
		self.databases = os.path.abspath(kwargs["dbdir"])
		self.plot = kwargs["plot"]
		samples = set()
		with smart_open(sample_f) as sam:
			for line in sam:
				if line.startswith("#"):
					continue
				sample_name = str(line.strip().split()[0])
				samples.add(sample_name)
				dep = os.path.join(depth_dir, "{0}/depthAnno/Flank/{0}.Rmdup.depth.tsv.gz".format(sample_name))
				self.SamplesDepth[sample_name] = dep
		if len(samples) <= 5:
			raise IOError('Sorry, Six or more samples are needed !!!')
		self.samples = sorted(samples)
		with smart_open(self.region) as bed:
			for line in bed:
				if line.startswith("#"):
					continue
				chrom_name = str(line.strip().split()[0])
				self.contigs.add(chrom_name)
		prepair_dir = FolderMaker(self.outdir)
		prepair_dir.create_subdirectory(self.samples)
		if not os.path.isdir(self.script):
			os.makedirs(self.script)

	def run(self):
		self.Step9_bedAnalysis()
		self.Step10_GCcorrect()
		self.Step11_WinCorrect()
		self.Step12_CNVCall()

	def Step9_bedAnalysis(self):
		jobs_tool = bedAnalysis(bedAnalysis_path=os.path.join(Rhea_Chip_Home, "tools", "cnv", "bedAnalysis.py"))
		s_out = os.path.join(self.script, "Step9.bedAnalysis.sh")
		bash = BASH("bedAnalysis", s_out, self.log)
		cd = jobs_tool.analysis(region=self.region, reference=self.reference, rmsk=self.rmsk,
		                        outdir=self.outdir, depthfile=",".join(self.SamplesDepth.values()))
		cd = bash.stat_command_formact(cd, "bed region analysis")
		bash.write(bash.bash_header_and_foot(cd, 'CNV bed region analysis'))
		bash.close()

	def Step10_GCcorrect(self):
		jobs_tool = GCcorrect(GCcorrect_path=os.path.join(Rhea_Chip_Home, "tools", "cnv", "GCcorrect.py"))
		wingc = os.path.join(self.outdir, "win.gc")
		posgc = os.path.join(self.outdir, "pos.gc")
		chromstat = os.path.join(self.outdir, "chrom.stat")
		for sample in self.samples:
			cd = jobs_tool.correct(sample=sample, input=self.SamplesDepth[sample],
			                       outdir=self.outdir, wingc=wingc,
			                       posgc=posgc, chromstat=chromstat)
			s_out = os.path.join(self.script, "Step10.GCcorrect.%s.sh" % sample)
			bash = BASH("%s.GCcorrect" % sample, s_out, self.log)
			cd = bash.stat_command_formact(cd, "Sample depth correct with GC")
			bash.write(bash.bash_header_and_foot(cd, 'Sample depth correct with GC'))
			bash.close()

	def Step11_WinCorrect(self):
		jobs_tool = WinCorrect(WinCorrect_path=os.path.join(Rhea_Chip_Home, "tools", "cnv", "WinCorrect.py"))
		for chrom in self.contigs:
			s_out = os.path.join(self.script, "Step11.WinCorrect.%s.sh" % chrom)
			bash = BASH("%s.WinCorrect" % chrom, s_out, self.log)
			cd = jobs_tool.analysis(region=self.region, chromstat=os.path.join(self.outdir, "chrom.stat"), chrom=chrom,
			                        indir=self.outdir, correct_win_len=self.winlen, correct_shift_len=self.siftlen)
			cd = bash.stat_command_formact(cd, "%s windows analysis and Correct" % chrom)
			bash.write(bash.bash_header_and_foot(cd, '%s windows analysis and Correct' % chrom))
			bash.close()

	def Step12_CNVCall(self):
		jobs_tool = batCNV(batCNV_path=os.path.join(Rhea_Chip_Home, "tools", "cnv", "batCNV.py"))
		for sample in self.samples:
			s_out = os.path.join(self.script, "Step12.CNVCall.%s.sh" % sample)
			bash = BASH("%s.CNVCall" % sample, s_out, self.log)
			cd = jobs_tool.analysis(indir=self.outdir, dbdir=self.databases, sample=sample, correct_win_len=self.winlen,
			                        correct_shift_len=self.siftlen, reference=self.reference, plot=self.plot)
			cd = bash.stat_command_formact(cd, "%s CNV calling" % sample)
			bash.write(bash.bash_header_and_foot(cd, '%s CNV calling' % sample))
			bash.close()


def main():
	usage = 'Usage: %prog [-h] [--version] --bed [region file] --samplelist [samplelist] [options]'
	description = cnv.__description__
	author = cnv.__author__
	version = cnv.__version__
	parser = OptionParser(usage=usage, version=version, description=description, epilog=author)
	expect = OptionGroup(parser, 'Expected arguments', 'Caution: These parameters are necessary for depthQC.')
	expect.add_option('-b', '--bed', metavar='FILE', dest='region', help='region file')
	expect.add_option('-l', '--samplelist', metavar='FILE', dest='samplelist', help='sample list file')
	parser.add_option_group(expect)
	optinal = OptionGroup(parser, 'Optional arguments', 'Caution: If you do not set these parameters in addition,'
	                                                    ' prog will select the default value.')
	optinal.add_option('-q', '--qc', metavar='DIR', dest='depths', default=None,
	                   help='QC dir, [ default: dir(sample.list)/QC ]')
	optinal.add_option('-o', '--outdir', dest='outdir', default=os.getcwd(),
	                   help='The output dir [ default: %s ]' % os.getcwd())
	optinal.add_option('-r', '--repeatdb', dest='repeatdb',
	                   default=os.path.join(Rhea_Chip_Home, "db/RepeatMasker/rmsk.bed.gz"),
	                   help='The RepeatMasker database dir [ default: %s ]' %
	                        os.path.join(Rhea_Chip_Home, "db/RepeatMasker/rmsk.bed.gz"))
	optinal.add_option('-d', '--database', dest='dbdir', default=os.path.join(Rhea_Chip_Home, "db"),
	                   help='The database dir [ default: %s ]' % os.path.join(Rhea_Chip_Home, "db"))
	optinal.add_option('-R', '--reference', dest='reference',
	                   default=os.path.join(Rhea_Chip_Home, "db", "aln_db", "hg19", "hg19_chM_male_mask.fa"),
	                   help='The human reference')
	optinal.add_option('-w', '--winlen', dest='winlen', default=30, help='The windowns length [ default: 30 ]')
	optinal.add_option('-s', '--siftlen', dest='siftlen', default=25, help='The windowns sift length [ default: 25 ]')
	optinal.add_option('-p', '--plot', help='Plot exon depth graph', dest='plot', default=False, action="store_true")
	parser.add_option_group(optinal)
	(options, args) = parser.parse_args()
	if not options.region or not options.samplelist:
		parser.print_help()
		return 'Expected arguments lost !!!'
	cnv_analysis = Pipeline(**options.__dict__)
	cnv_analysis.run()


if __name__ == '__main__':
	sys.exit(main())
