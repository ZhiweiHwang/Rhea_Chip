#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
from collections import defaultdict


class bedAnalysis(object):
	def __init__(self, **kwargs):
		kwargs = defaultdict(lambda: None, kwargs)
		tool_path = kwargs["bedAnalysis_path"]
		self.tools = os.path.abspath(tool_path) if tool_path else "bedAnalysis.py"

	def analysis(self, **kwargs):
		kwargs = defaultdict(lambda: None, kwargs)
		bed = kwargs["region"]
		reference = kwargs["reference"]
		depths = kwargs["depthfile"]
		options = ' -b %s -r %s -f %s' % (bed, reference, depths)
		if kwargs["rmsk"]:
			options += " -d %s" % kwargs["rmsk"]
		if kwargs["outdir"]:
			options += " -o %s" % kwargs["outdir"]
		if kwargs["winlen"]:
			options += " -w %i" % int(kwargs["winlen"])
		if kwargs["siftlen"]:
			options += " -s %i" % int(kwargs["siftlen"])
		return "{0} {1}".format(self.tools, options)
