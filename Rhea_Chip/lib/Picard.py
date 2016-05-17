#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
from JavaTools import JavaTool
from Samtools import Samtools


class Picard(JavaTool):
	def __init__(self, jar_path="", javatmp=os.getcwd(), max_memory=None, **kwargs):
		if not jar_path:
			raise IOError("Fail to load picard tools~")
		JavaTool.__init__(self, os.path.abspath(jar_path), tmp_path=os.path.abspath(javatmp), max_memory=max_memory)

	def SortSam(self, in_bam, out_bam=None, sort_order='coordinate', **kwargs):
		jar = 'SortSam.jar'
		in_bam = os.path.abspath(in_bam)
		in_prefix, _ = os.path.splitext(in_bam)
		out_bam = os.path.abspath(out_bam) if out_bam else in_prefix + ".sort.bam"
		return "{0} INPUT={1} OUTPUT={2} VALIDATION_STRINGENCY=SILENT SORT_ORDER={3}".format(
				self.options + '/' + jar, in_bam, out_bam, sort_order)

	def MarkDuplicates(self, in_bam, out_bam=None, max_file_handles=8000, metrics_file=None, **kwargs):
		jar = 'MarkDuplicates.jar'
		in_bam = os.path.abspath(in_bam)
		in_prefix, _ = os.path.splitext(in_bam)
		out_bam = os.path.abspath(out_bam) if out_bam else in_prefix + ".sort.dup.bam"
		metrics_file = os.path.abspath(metrics_file) if metrics_file else in_prefix + ".dup.metrics"
		max_file_handles = "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=%i" % int(max_file_handles) if max_file_handles else ""
		return "{0} INPUT={1} OUTPUT={2} VALIDATION_STRINGENCY=SILENT METRICS_FILE={3} {4}".format(
				self.options + '/' + jar, in_bam, out_bam, metrics_file, max_file_handles)

	def FixMateInformation(self, in_bam, out_bam=None, **kwargs):
		jar = 'FixMateInformation.jar'
		in_bam = os.path.abspath(in_bam)
		in_prefix, _ = os.path.splitext(in_bam)
		out_bam = os.path.abspath(out_bam) if out_bam else in_prefix + ".fix_mate_infor.bam"
		return "{0} INPUT={1} OUTPUT={2} VALIDATION_STRINGENCY=SILENT".format(self.options + '/' + jar, in_bam, out_bam)

	def MergeSamFiles(self, input, output, **kwargs):
		jar = "MergeSamFiles.jar"
		bamlist = input.split()
		if len(bamlist) < 1:
			raise IOError("no bam file found ~")
		output = os.path.abspath(output)
		if len(bamlist) == 1:
			samtools = Samtools()
			if input.endswith('.sam') and output.endswith('.bam'):
				return samtools.view(input_file=input ,output_file=output, output_as_sam=True)
			elif input.endswith('.bam') and output.endswith('.sam'):
				return samtools.view(input_file=input ,output_file=output, include_header_in_output=True)
			else:
				return "cp {0} {1}".format(input, output)
		in_bam = " ".join(['INPUT=%s' % os.path.abspath(bam) for bam in bamlist])
		return "{0} {1} OUTPUT={2} VALIDATION_STRINGENCY=SILENT".format(self.options + '/' + jar, in_bam, output)
