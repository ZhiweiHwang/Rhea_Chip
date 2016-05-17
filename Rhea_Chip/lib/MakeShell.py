#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！
import os
import sys
from smart_open import smart_open


Rhea_Chip_Home = os.getenv('Rhea_Chip_Home')


class BASH(object):
	def __init__(self, keyword, output=None, logfile=None):
		if output:
			self.output = smart_open(os.path.abspath(output), 'w')
		else:
			self.output = sys.stdout
		self.keyword = keyword.strip()
		self.log = ">> %s" % os.path.abspath(logfile) if logfile else ""

	def stat_command_formact(self, messages, keyword, level=1):
		cd = [messages, 'if [ $? -ne 0 ]; then', '\techo "{0} {1} failed." {2}'.format(self.keyword, keyword, self.log),
		      '\texit 100', 'fi\n']
		levels = "\n" + "\t" * (level - 1)
		return "\t" * (level - 1) + levels.join(cd)

	def bash_header_and_foot(self, meassagelist, step):
		if not len(meassagelist):
			return list()
		commandlist = ['#!/bin/bash\n',
		               'program_start=`date "+%s"`\n',
		               'export PATH=%s/tools:$PATH\n' % Rhea_Chip_Home,
		               '\nprogram_end=`date "+%s"`\n',
		               'run_time=`expr "$program_end" - "$program_start"`\n',
		               'echo - | awk -v S=$run_time \'{printf "%s %s spends %%02d:%%02d:%%02d\\n",'
		               ' S/(60*60),S%%(60*60)/60,S%%60}\' %s\nexit 0\n' % (self.keyword, step, self.log)]
		if type(meassagelist) is list:
			commandlist[3:3] = meassagelist
		elif type(meassagelist) is str:
			commandlist.insert(3, meassagelist)
		elif type(meassagelist) is tuple:
			commandlist[3:3] = list(meassagelist)
		else:
			return list()
		return commandlist

	def write(self, lines, **kwargs):
		if type(lines) is str:
			self.output.write(lines + '\n')
		else:
			self.output.writelines(lines)

	def close(self):
		self.output.close()