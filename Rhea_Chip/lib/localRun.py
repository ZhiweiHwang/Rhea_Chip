#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！
import multiprocessing
from subprocess import call


class Local(object):
	def __init__(self, count=None):
		c_new = multiprocessing.cpu_count()
		count = int(count) if count else c_new
		self.pool = multiprocessing.Pool(processes=min(count, c_new))

	def run(self, script, err=None, log=None, **kwargs):
		if not err:
			err = script + '.err'
		if not log:
			err = script + '.out'
		script = "sh %s 2>%s >%s" % (script, err, log)
		self.pool.apply_async(work, (script,))

	def __del__(self):
		self.pool.close()
		self.pool.join()

	def close(self):
		self.__del__()

def work(cmd):
	return call(cmd.split(), shell=False)