#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
import re
import Rhea_Chip
import sys
from itertools import groupby
from collections import defaultdict
from glob import glob
from optparse import OptionParser
from Rhea_Chip.lib import Qsub, Local
from Rhea_Chip.lib.ParserConfig import ParserConfig


def scripts_phrase(scrips_d, samples=set()):
	scrips_dict = dict()
	for i in filter(lambda s: os.path.isdir(os.path.join(scrips_d, s)), os.listdir(scrips_d)):
		if len(samples) and i not in samples:
			continue
		scrips_dict[i] = os.path.abspath(os.path.join(scrips_d, i))
	return scrips_dict


class Run(object):
	def __init__(self, model="local", **kwargs):
		kwargs = defaultdict(lambda: None, kwargs)
		self.model = model
		if self.model == "sge":
			self.run = Qsub(queue_parameter=kwargs['queue_parameter'])
		elif self.model == "local":
			self.run = Local(count=kwargs['count'])
		else:
			raise "model %s not support !" % self.model

	def start(self, scripts, hold_jid=None):
		if self.model == "sge":
			orders = dict()
			for (s, m) in sorted(scripts, key=lambda x: self.order_s(x[0])):
				if not os.path.isfile(s):
					continue
				err = re.sub("\.sh$", ".err", s)
				log = re.sub("\.sh$", ".out", s)
				num = self.order_s(s)
				try:
					hold = ",".join([str(i) for i, j in orders.iteritems() if int(j) < num]) if len(orders) else str()
					hold += ",%s" % hold_jid if hold_jid else ""
					ids, jobs = self.run.run(s, mem=m, err=err, log=log, hold=hold)
					orders[ids] = num
					print "ID [ {0} ] : Task [ {1} ] has been successfully delivered !".format(ids, jobs)
				except IOError as e:
					self.run.close()
					raise IOError(e)
			return orders
		elif self.model == "local":
			orders = dict()
			for k, g in groupby(scripts, key=lambda x: self.order_s(x)):
				orders[k] = list(g)
			for s, m in sorted(orders.iteritems(), key=lambda x: x[0]):
				for j in m:
					err = re.sub("\.sh$", ".err", j)
					log = re.sub("\.sh$", ".out", j)
					self.run.run(j, err, log)
				self.run.close()
			return orders
		else:
			pass

	@staticmethod
	def order_s(name):
		name = os.path.basename(name)
		order = re.compile("step([0-9]+)\..*\.sh$", re.I).match(name)
		if order:
			return int(order.group(1))
		else:
			return 0


def main():
	usage = 'Usage: %prog [-h] [--version] --samplelist [sample list] --config [config file] --dir [scripts dir]'
	description = Rhea_Chip.__description__
	author = Rhea_Chip.__author__
	version = Rhea_Chip.__version__
	parser = OptionParser(usage=usage, version=version, description=description, epilog=author)
	parser.add_option('-s', '--samplelist', metavar='FILE', dest='samplelist', default=None,
	                  help=r'''The sample list file
	                        Column 1:    sample name
	                        Column 2:    family name
	                        Column 3:    library
	                        Column 4:    FQ file [formact:
	                                                 lane1/fq1:lane1/adapter1;lane1/fq2:lane1/adapter2|
	                                                 lane2/fq1:lane2/adapter1;lane2/fq2:lane2/adapter2|...]
	                        Column 5:    case/control''')
	parser.add_option('-c', '--config', metavar='FILE', dest='config',
	                  help='The config file, contain all specific settings')
	parser.add_option('-d', '--dir', dest='dirs', default=None,
	                  help='The scripts dir, contain all of the scripts of all samples')
	parser.add_option('-C', '--CNV', dest='cnv_dirs', default=None,
	                  help='The cnv scripts dir, contain all of the cnv calling scripts of all samples')
	options, args = parser.parse_args()
	if not (options.samplelist or options.dirs) or not options.config:
		parser.print_help()
		return 'Expected arguments lost !!!'
	samples = options.samplelist
	sample_set = set()
	if samples is not None and os.path.isfile(samples):
		with open(samples) as f_in:
			for line in f_in:
				sample = line.strip().split()[0]
				sample_set.add(sample)
	scrips_d = os.path.abspath(options.dirs) if options.dirs else os.path.abspath(os.path.join(samples, "../../script"))
	sample_dict = scripts_phrase(scrips_d, sample_set)
	if not len(sample_set):
		sample_set = set(sample_dict.keys())
	baseconfig = ParserConfig(options.config)
	cnv_dir = os.path.abspath(options.cnv_dirs) if options.cnv_dirs else None
	job = baseconfig.handle.get("Analysis Flow", "analysis_option") or baseconfig.handle.get('Analysis Flow',
	                                                                                         'Rhea_Flow_analysis')
	jobs = job.split(",")
	run_model = baseconfig.handle.get("Analysis Flow", "tasks_management_option") or "local"
	queue_parameter = baseconfig.handle.get("Qsub Parameter", "queue_parameter")
	run = Run(model=run_model, queue_parameter=queue_parameter, count=len(sample_set))
	all_jobs = defaultdict(set)

	for sample in sample_set:
		if sample not in sample_dict:
			raise IOError("The sample %s was not found in the scripts' directory !" % sample)
		job_set = set()
		scripts = glob(os.path.join(sample_dict[sample], "Step*.sh"))
		for script in scripts:
			job_order, job_name = os.path.basename(script).split(".")[:2]
			if job_name in jobs and job_name != "result":
				jobs_mem = baseconfig.handle.get("Qsub Parameter", ".".join([job_order, job_name])) or "2G"
				job_set.add((script, jobs_mem))
		orders = run.start(job_set)
		for i, j in orders.iteritems():
			all_jobs[j].add(i)

	if 'cnv' in jobs and cnv_dir is not None:
		job_set = set()
		hold_jobs = list()
		start_order = sys.maxint
		scripts = glob(os.path.join(cnv_dir, "Step*.sh"))
		for script in scripts:
			job_order, job_name = os.path.basename(script).split(".")[:2]
			start_order = min(run.order_s(script), start_order)
			jobs_mem = baseconfig.handle.get("Qsub Parameter", ".".join([job_order, job_name])) or "2G"
			job_set.add((script, jobs_mem))
		for i, j in all_jobs.iteritems():
			if i < start_order:
				hold_jobs.extend(map(str, list(j)))
		orders = run.start(job_set, hold_jid=",".join(hold_jobs))
		for i, j in orders.iteritems():
			all_jobs[j].add(i)

	if "result" in jobs:
		scripts = glob(os.path.join(scrips_d, "*", "Step*.result.sh"))
		start_order = run.order_s(scripts[0])
		job_name = ".".join(os.path.basename(scripts[0]).split(".")[:2])
		jobs_mem = baseconfig.handle.get("Qsub Parameter", job_name) or "2G"
		job_set = zip(scripts, [jobs_mem, ] * len(scripts))
		hold_jobs = list()
		for i, j in all_jobs.iteritems():
			if i < start_order:
				hold_jobs.extend(map(str, list(j)))
		_ = run.start(job_set, hold_jid=",".join(hold_jobs))


if __name__ == '__main__':
	sys.exit(main())
