#! /usr/bin/env python
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
from collections import defaultdict


class GCcorrect(object):
	def __init__(self, **kwargs):
		kwargs = defaultdict(lambda: None, kwargs)
		tool_path = kwargs["GCcorrect_path"]
		self.tools = os.path.abspath(tool_path) if tool_path else "GCcorrect.py"

	def correct(self, **kwargs):
		kwargs = defaultdict(lambda: None, kwargs)
		wingc = os.path.abspath(kwargs["wingc"])
		posgc = os.path.abspath(kwargs["posgc"])
		inputs = os.path.abspath(kwargs["input"])
		chromstat = os.path.abspath(kwargs["chromstat"])
		options = ' -w %s -p %s -i %s -c %s' % (wingc, posgc, inputs, chromstat)
		if kwargs["outdir"]:
			options += " -o %s" % kwargs["outdir"]
		if kwargs["sample"]:
			options += " -s %s" % str(kwargs["sample"])
		return "{0} {1}".format(self.tools, options)



