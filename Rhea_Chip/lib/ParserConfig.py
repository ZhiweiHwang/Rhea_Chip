#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
from sys import argv
from ConfigParser import ConfigParser


class ParserConfig(object):
	def __init__(self, configfile):
		self.config = os.path.abspath(configfile)
		self.handle = ConfigParser()
		self.handle.read(self.config)

	def config_check(self):
		try:
			flows = set(self.handle.get("Analysis Flow", 'Rhea_flow_analysis').split(","))
			analysis = set(self.handle.get("Analysis Flow", 'analysis_option').split(","))
			new_job = analysis.difference(flows)
			if len(new_job):
				return "Sorry, no support Operations [ %s ] yet !" % ",".join(new_job)
			male_ref = self.handle.get("Reference databases", 'male_ref')
			female_ref = self.handle.get("Reference databases", 'female_ref')
			if not (os.path.isfile(male_ref) and os.path.isfile(male_ref + '.fai') and os.path.isfile(
					female_ref) and os.path.isfile(female_ref + '.fai')):
				return "Sorry, Reference databases must be set, and all of them should be cread index ~"
			target_region = self.handle.get("Personal Analysis", "target_region")
			run_region = self.handle.get("Personal Analysis", "run_region")
			if not os.path.isfile(target_region) and not os.path.isfile(run_region):
				return "Sorry, you must define more than one bed file ~"
			if not os.path.isfile(target_region):
				self.handle.set("Personel Analysis", "target_region", run_region)
			if not os.path.isfile(run_region):
				self.handle.set("Personel Analysis", "run_region", target_region)
		except Exception:
			return "Config file damaged !!!"
		return "Success"

	def config_save(self, config_out):
		config_out = os.path.abspath(config_out)
		if config_out == self.config:
			return
		f_out = open(config_out, "w")
		self.handle.write(f_out)
		f_out.close()
		return


if __name__ == '__main__':
	pro = argv[0]
	if len(argv) < 2:
		raise IOError("python %s <f_in> [f_out]" % pro)
	f_in = argv[1]
	cf = ParserConfig(f_in)
	print "Config file check result : %s" % cf.config_check()
	if len(argv) > 2:
		f_out = argv[2]
		cf.config_save(f_out)
