#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os


class Samtools(object):
	def __init__(self, samtools_path=None, max_threads=0, **kwargs):
		self.tools = os.path.abspath(samtools_path) if samtools_path else "samtools"
		self.max_threads = int(max_threads)

	def index(self, in_bam, out=None, **kwargs):
		in_bam = os.path.abspath(in_bam)
		out = os.path.abspath(out) if out else in_bam + '.bai'
		return "{0} index {1} {2}".format(self.tools, in_bam, out)

	def rmdup(self, input_bam, output_bam, remove_dup_for_se_reads=False, treat_both_pe_and_se_reads=False, **kwargs):
		options = " -s" if remove_dup_for_se_reads else ""
		options += " -S" if treat_both_pe_and_se_reads else ""
		options += " %s" % os.path.abspath(input_bam)
		options += " %s" % os.path.abspath(output_bam)
		return "{0} rmdup {1}".format(self.tools, options)

	def view(self, input_file=None, output_file=None, include_header_in_output=False, output_as_sam=False,
	         white_flag_value=None, black_flag_value=None, bed_file_with_regions_to_output=None,
	         auto_detected=True, **kwargs):
		options = " -h" if include_header_in_output else ""
		options += "" if output_as_sam else " -b"
		options += " -S" if auto_detected else ""
		options += " -f %i" % int(white_flag_value) if white_flag_value is not None else ""
		options += " -F %i" % int(black_flag_value) if black_flag_value is not None else ""
		options += " -L %s" % bed_file_with_regions_to_output if bed_file_with_regions_to_output else ""
		options += " %s" % os.path.abspath(input_file) if input_file else " -"
		options += " -o %s" % os.path.abspath(output_file) if output_file else "|"
		return "{0} view {1}".format(self.tools, options)

	def sort(self, input_bam, output_bam, temp_file_prefix=None, **kwargs):
		if temp_file_prefix is None:
			temp_file_prefix = output_bam + '.tmp.bam'
		options = " -@ %i" % self.max_threads if self.max_threads else ""
		options += " -o %s" % os.path.abspath(output_bam)
		options += " -T %s" % os.path.abspath(temp_file_prefix)
		options += " %s" % os.path.abspath(input_bam)
		return "{0} sort {1}".format(self.tools, options)

	def faidx(self, fasta_file, **kwargs):
		return "{0} faidx {1}".format(self.tools, os.path.abspath(fasta_file))
