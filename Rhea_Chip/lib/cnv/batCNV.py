#! /usr/bin/env python
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！
import os
from collections import defaultdict


class batCNV(object):
	def __init__(self, **kwargs):
		kwargs = defaultdict(lambda: None, kwargs)
		tool_path = kwargs["batCNV_path"]
		self.tools = os.path.abspath(tool_path) if tool_path else "batCNV.py"

	def analysis(self, **kwargs):
		kwargs = defaultdict(lambda: None, kwargs)
		indir = kwargs["indir"]
		dbdir = kwargs["dbdir"]
		options = ' -i %s -d %s' % (indir, dbdir)
		if kwargs["correct_win_len"]:
			options += " --correct_win_len %i" % int(kwargs["correct_win_len"])
		if kwargs["correct_shift_len"]:
			options += " --correct_shift_len %i" % int(kwargs["correct_shift_len"])
		if kwargs["reference"]:
			options += " -r %s" % os.path.abspath(kwargs["reference"])
		if kwargs["chrom"]:
			options += " -c %s" % kwargs["chrom"]
		if kwargs["sample"]:
			options += " -s %s" % str(kwargs["sample"])
		if kwargs["plot"]:
			plot = True if str(kwargs["plot"]) == 'True' else False
			if plot:
				options += " -p"
		return "{0} {1}".format(self.tools, options)
