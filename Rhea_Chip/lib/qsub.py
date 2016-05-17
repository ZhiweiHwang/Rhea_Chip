#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！
import os
import re
from subprocess import Popen, PIPE


class Qsub(object):
	def __init__(self, tools=None, queue_parameter=None, **kwargs):
		self.qsub = os.path.abspath(tools) if tools else "qsub"
		self.queue_parameter = queue_parameter if queue_parameter else ""
		self.qdel = os.path.join(os.path.dirname(self.qsub), "qdel") if tools else "qdel"
		self.qmod = os.path.join(os.path.dirname(self.qsub), "qmod") if tools else "qmod"
		self.jobs = dict()

	def run(self, script, mem="2g", err=None, log=None, hold=None, **kwargs):
		cmd = "{0} -cwd -l vf={1} {2} ".format(self.qsub, mem, self.queue_parameter)
		if err:
			cmd += "-e %s " % os.path.abspath(err)
		if log:
			cmd += "-o %s " % os.path.abspath(log)
		if hold:
			cmd += "-hold_jid %s " % hold
		cmd += os.path.abspath(script)
		out, err = Popen(cmd.split(), stdout=PIPE, stderr=PIPE).communicate()
		try:
			job_id = re.compile("^Your job ([0-9]+) \(\"([^\)]*)\"\) has been submitted$").search(out.strip()).groups()
			ids = int(job_id[0])
			scripts = os.path.basename(job_id[1])
			self.jobs[scripts] = str(ids)
			return ids, scripts
		except:
			raise IOError("%s\nExecuted command:\n%s\n" % (err, cmd))

	def close(self, jobs=None, **kwargs):
		jobs = str(jobs) if jobs else " ".join(self.jobs.values())
		if jobs:
			cmd = [self.qdel, jobs]
			Popen(cmd, stdout=PIPE, stderr=PIPE).communicate()
		else:
			return

	def mod(self, jobs=None, **kwargs):
		jobs = str(jobs) if jobs else " ".join(self.jobs.values())
		if jobs:
			cmd = [self.qmod, '-cj', jobs]
			Popen(cmd, stdout=PIPE, stderr=PIPE).communicate()
		else:
			return

	def __iter__(self):
		return sorted(self.jobs.iteritems(), key=lambda s: int(s[1]))
