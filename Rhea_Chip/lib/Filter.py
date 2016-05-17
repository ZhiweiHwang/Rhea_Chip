#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os


class SOAPnuke(object):
	def __init__(self, soapnuke_path=None, **kwargs):
		self.tools = os.path.abspath(soapnuke_path) if soapnuke_path else "SOAPnuke"

	def filter(self, adapter1=None, adapter2=None, fq1=None, fq2=None, outdir=os.getcwd(), clean_fq1=None,
	           clean_fq2=None, clean_fq_type='sanger', qual_threshold=10, ave_qual_filter=0, qual_filter_rate=0.5,
	           n_rate=0.1, qual_sys='illumina', seq_type='new', poly_a_filter=0, rmdup=0, **kwargs):
		outdir = os.path.abspath(outdir)
		options = "-o %s -l %i -q %f -n %f" % (outdir, int(qual_threshold), float(qual_filter_rate), float(n_rate))
		if fq1 and clean_fq1:
			adapter1 = "-f %s" % os.path.abspath(adapter1) if adapter1 else ""
			fq1 = os.path.abspath(fq1)
			clean_fq1 = os.path.basename(clean_fq1)
			options += " %s -1 %s -C %s" % (adapter1, fq1, clean_fq1)
		if fq2 and clean_fq2:
			adapter2 = "-r %s" % os.path.abspath(adapter2) if adapter2 else ""
			fq2 = os.path.abspath(fq2)
			clean_fq2 = os.path.basename(clean_fq2)
			options += " %s -2 %s -D %s" % (adapter2, fq2, clean_fq2)
		if float(ave_qual_filter) > 0:
			options += " -m %f" % float(ave_qual_filter)
		options += " -G" if clean_fq_type == 'sanger' else ""
		options += " -Q 2" if qual_sys == 'sanger' else " -Q 1 "
		options += " -5 1" if seq_type == 'new' else ""
		options += " -p %f" % float(poly_a_filter) if float(poly_a_filter) > 0 else " -p 0 "
		if int(rmdup):
			options += " -d"
		return "{0} filter {1}".format(self.tools, options)
