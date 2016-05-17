#! /usr/bin/env python
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
from collections import defaultdict


class WinCorrect(object):
	def __init__(self, **kwargs):
		kwargs = defaultdict(lambda: None, kwargs)
		tool_path = kwargs["WinCorrect_path"]
		self.tools = os.path.abspath(tool_path) if tool_path else "WinCorrect.py"

	def analysis(self, **kwargs):
		kwargs = defaultdict(lambda: None, kwargs)
		bed = kwargs["region"]
		chromstat = kwargs["chromstat"]
		indir = kwargs["indir"]
		options = ' -b %s -c %s -i %s' % (bed, chromstat, indir)
		if kwargs["low_dep_cut"]:
			options += " --low_dep_cut %.3f" % float(kwargs["low_dep_cut"])
		if kwargs["correct_win_len"]:
			options += " --correct_win_len %i" % int(kwargs["correct_win_len"])
		if kwargs["correct_shift_len"]:
			options += " --correct_shift_len %i" % int(kwargs["correct_shift_len"])
		if kwargs["chrom"]:
			options += " -C %s " % str(kwargs["chrom"])
		return "{0} {1}".format(self.tools, options)