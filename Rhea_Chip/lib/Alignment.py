#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
from collections import defaultdict

class BWA(object):
	def __init__(self, **kwargs):
		kwargs = defaultdict(lambda :None, kwargs)
		bwa_path = kwargs["bwa_path"]
		max_threads = kwargs["max_threads"]
		self.tools = os.path.abspath(bwa_path) if bwa_path else "bwa"
		self.max_threads = int(max_threads) if max_threads else 4

	def index(self, reference, **kwargs):
		kwargs = defaultdict(lambda :None, kwargs)
		reference = os.path.abspath(reference)
		prefix = kwargs["prefix"]
		index_prefix = "-p %s" % os.path.abspath(prefix) if prefix else ""
		return "{0} index {1} {2}".format(self.tools, index_prefix, reference)

	def mem(self, reads, reference, other_options="", **kwargs):
		kwargs = defaultdict(lambda :None, kwargs)
		readgroup = kwargs["readgroup"]
		output = kwargs["output"]
		reference = os.path.abspath(reference)
		if readgroup:
			other_options += " -R \"%s\"" % readgroup
		output = "> %s" % os.path.abspath(output) if output else "|"
		return "{0} mem -t {1} {2} {3} {4} {5}".format(self.tools, self.max_threads,
		                                               other_options, reference, reads, output)

	def aln(self, reads, reference, out_sai, other_options="", **kwargs):
		output = "-f %s" % os.path.abspath(out_sai)
		return "{0} aln -t {1} {2} {3} {4} {5}".format(self.tools, self.max_threads,
		                                               other_options, output, reference, reads)

	def samse(self, in_sai, in_reads, reference, **kwargs):
		kwargs = defaultdict(lambda :None, kwargs)
		readgroup = kwargs["readgroup"]
		output = kwargs["output"]
		max_occ = kwargs["max_occ"]
		other_options = "-n %i" % int(max_occ) if max_occ else ""
		output = "-f %s" % os.path.abspath(output) if output else "|"
		reference = os.path.abspath(reference)
		if readgroup:
			other_options += "-r \"%s\"" % readgroup
		return "{0} samse {1} {2} {3} {4} {5}".format(self.tools, other_options, reference, in_sai, in_reads, output)

	def sampe(self, in_sai, in_reads, reference, max_insert_size=500, max_occ=100000, max_hit_pair=3,
	          max_hit_discordant=10, **kwargs):
		kwargs = defaultdict(lambda :None, kwargs)
		readgroup = kwargs["readgroup"]
		output = kwargs["output"]
		options = "-a %i -o %i -n %i -N %i" % (max_insert_size, max_occ, max_hit_pair, max_hit_discordant)
		output = "-f %s" % os.path.abspath(output) if output else "|"
		reference = os.path.abspath(reference)
		if readgroup:
			options += "-r \"%s\"" % readgroup
		return "{0} sampe {1} {2} {3} {4} {5}".format(self.tools, options, reference, in_sai, in_reads, output)
