#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
import pysam
import tempfile
from Rhea_Chip.lib.CreateFolder import mksymlink


temp_directory_name = tempfile.mkdtemp()


class CreatIndex(object):
	def __init__(self, files, indexdir=temp_directory_name):
		outdir, self.filename = os.path.split(files)
		self.indexdir = os.path.abspath(indexdir)
		self.files = os.path.abspath(files)
		self.outdir = os.path.abspath(outdir)
		self.prefix, self.ext = os.path.splitext(files)

	def check_index(self, preset=None, seq_col=None, start_col=None, end_col=None):
		if self.ext == 'bam' or preset == 'bam':
			file_index = self.files + '.bai'
			if not os.path.isfile(file_index):
				file_index = self.prefix + '.bai'
			if os.path.isfile(file_index) and os.path.getmtime(file_index) >= os.path.getmtime(self.files):
				pass
			else:
				self.files = self.make_index()
				_ = pysam.index(self.files, force=True)
		elif self.ext == 'fa' or preset == 'fasta' or self.ext == 'fasta':
			file_index = self.files + '.fai'
			if os.path.isfile(file_index) and os.path.getmtime(file_index) >= os.path.getmtime(self.files):
				pass
			else:
				self.files = self.make_index()
				_ = pysam.FastaFile(self.files)
		elif self.ext == 'vcf' or preset == 'vcf' or self.files.endswith('vcf.gz'):
			file_index = self.files + '.tbi'
			if os.path.isfile(file_index) and os.path.getmtime(file_index) >= os.path.getmtime(self.files):
				pass
			else:
				self.files = self.make_index()
				self.files = pysam.tabix_index(self.files, preset='vcf', force=True)
		else:
			file_index = self.files + '.tbi'
			if os.path.isfile(file_index) and os.path.getmtime(file_index) >= os.path.getmtime(self.files):
				pass
			else:
				self.files = self.make_index()
				try:
					self.files = pysam.tabix_index(self.files, preset=preset, seq_col=seq_col,
					                               start_col=start_col, end_col=end_col, force=True)
				except Exception as err:
					print Exception(err)
		return self.files

	def make_index(self):
		if not os.access(self.outdir, os.W_OK):
			mksymlink(self.files, os.path.join(self.indexdir, self.filename))
			files = os.path.join(self.indexdir, self.filename)
		else:
			files = self.files
		return files
