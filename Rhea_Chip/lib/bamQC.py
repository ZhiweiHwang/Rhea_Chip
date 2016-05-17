#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
from collections import defaultdict


class DepthQC(object):
	def __init__(self, **kwargs):
		path = kwargs["path"] if "path" in kwargs else None
		self.tools = os.path.abspath(path) if path else "depthQC.py"

	def run(self, in_bam, target_bed, reference, outdir=os.getcwd(), plot=True, max_threads=0, **kwargs):
		kwargs = defaultdict(lambda :None, kwargs)
		flankbed = kwargs["flankbed"]
		trans = kwargs["trans"]
		genes = kwargs["genes"]
		options = " -t %i" % int(max_threads) if int(max_threads) else ""
		options += " -b %s -o %s --target %s -r %s" % (in_bam, outdir, target_bed, reference)
		if flankbed:
			options += " --flank %s" % flankbed
		if trans:
			options += " --trans %s" % trans
		elif genes:
			options += " --genes %s" % genes
		else:
			raise KeyError("Genes and Trans must set one ~")
		if plot:
			options += " --plot "
		return self.tools + options


