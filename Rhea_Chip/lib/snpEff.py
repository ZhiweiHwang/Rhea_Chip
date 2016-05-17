#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
from JavaTools import JavaTool, defaultdict


class snpEff(JavaTool):
	def __init__(self, jar="snpEff.jar", javatmp=os.getcwd(), max_memory=None):
		if not jar:
			raise IOError("Fail to load snpEff tools~")
		self.jar = os.path.abspath(jar)
		JavaTool.__init__(self, self.jar, tmp_path=os.path.abspath(javatmp), max_memory=max_memory)

	def VariantAnnotator(self, vcf, **kwargs):
		vcf = os.path.abspath(vcf)
		kwargs = defaultdict(lambda :None, kwargs)
		mode = kwargs["mode"] if kwargs["mode"] else "hg19"
		output = "> %s" % os.path.abspath(kwargs["output"]) if kwargs["output"] else "|"
		config = os.path.abspath(kwargs["config"]) if kwargs["config"] else os.path.join(os.path.dirname(self.jar), "snpEff.config")
		if not os.path.isfile(config):
			raise IOError("Fail to load snpEff config file~")
		options = ' -c %s ' % config
		if kwargs["formats"]:
			options += "-o %s" % kwargs["formats"]
		if kwargs["canonical"]:
			options += "-canon "
		if kwargs["interval"] and os.path.isfile(kwargs["interval"]):
			options += "-interval %s " % os.path.abspath(kwargs["interval"])
		if kwargs["validated"]:
			options += "-strict "
		if kwargs["transfile"] and os.path.isfile(kwargs["transfile"]):
			options += "-onlyTr %s " % os.path.abspath(kwargs["transfile"])
		if kwargs["stats"]:
			stats = os.path.abspath(kwargs["stats"])
		else:
			stats_dir = os.path.dirname(os.path.abspath(kwargs["output"])) if kwargs["output"] else os.path.dirname(vcf)
			stats = os.path.join(stats_dir, "snpEff.html")
		options += "-s %s " % stats
		return "{0} {1} {2} {3} {4}".format(self.options, options, mode, vcf, output)
