#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
import time
import shutil
import re
from sys import argv


class FolderMaker(object):
	"""
	a simple class use for makedir
	"""

	def __init__(self, rootdir):
		self.rootdir = self.folder_confirm(os.path.abspath(rootdir))

	@staticmethod
	def folder_confirm(realdir):
		drive = realdir.split(os.sep)
		temp_dir = '/'
		while drive:
			temp_dir = os.path.join(temp_dir, drive.pop(0))
			try:
				os.makedirs(temp_dir)
			except OSError:
				if os.path.isdir(temp_dir):
					pass
				else:
					raise
		return realdir

	def create_subdirectory(self, sublist):
		temptime = time.strftime('%Y.%m.%d', time.localtime(time.time()))
		for dirpath in sublist:
			if not len(dirpath):
				continue
			newpath = os.path.join(self.rootdir, dirpath)
			try:
				if os.path.isfile(newpath):
					os.rename(newpath, str(newpath) + temptime)
				os.makedirs(newpath)
			except OSError:
				if os.path.isdir(newpath):
					pass
				else:
					raise
		return True


def mksymlink(oldpath, linkpath):
	"""
	Establish a soft link for the target file
	"""
	if not oldpath or not len(oldpath) or not os.path.exists(oldpath):
		return
	suffix = re.compile('(\\.[a-z]+)$', re.I).search(oldpath)
	if suffix and not linkpath.endswith(suffix.group()):
		linkpath += suffix.group()
	if os.path.isdir(linkpath):
		shutil.rmtree(linkpath)
	elif os.path.isfile(linkpath):
		os.remove(linkpath)
	os.symlink(oldpath, linkpath)


if __name__ == '__main__':
	if len(argv) < 3:
		raise IOError("python %s <root_dir> <sub_dir_list>" % argv[0])
	prog, root_dir, sub_dir = argv
	cf = FolderMaker(os.path.abspath(root_dir))
	cf.create_subdirectory(sub_dir)
