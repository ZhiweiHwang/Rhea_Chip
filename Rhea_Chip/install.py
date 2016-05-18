#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# The MIT License (MIT)
# Copyright (c) 2015 Joey Hwong
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
# associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions: The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
Install RheaChip pipeline.
"""

import os
import re
import shutil
import sys
import logging
from subprocess import call, PIPE, Popen
import time
import tarfile
from distutils.version import LooseVersion
from optparse import OptionParser
from collections import defaultdict

Rhea_Chip_Home = os.getenv('Rhea_Chip_Home') or os.path.dirname(os.path.abspath(__file__))
logger = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s: %(levelname)s: %(message)s', level=logging.INFO)
logger.info('Running: python %s' % ' '.join(sys.argv))


def mksymlink(oldpath, linkpath):
	if not oldpath or not len(oldpath) or not os.path.exists(oldpath):
		return False
	suffix = re.compile('(\\.[a-z]+)$', re.I).search(oldpath)
	if suffix and not linkpath.endswith(suffix.group()):
		linkpath += suffix.group()
	if os.path.isdir(linkpath):
		shutil.rmtree(linkpath)
	elif os.path.exists(linkpath):
		os.remove(linkpath)
	os.symlink(oldpath, linkpath)


class Install(object):
	def __init__(self, **kwargs):
		kwargs = defaultdict(lambda: None, kwargs)
		self._download_url = r'http://115.238.230.165:8010/get'
		self.prefix = os.path.abspath(kwargs['prefix']) if kwargs['prefix'] else Rhea_Chip_Home
		self.nodata = True if str(kwargs['nodata']) == "True" else False
		self.nosoft = True if str(kwargs['nosoft']) == "True" else False
		self.max_retries = 5
		self.files = defaultdict(dict)
		self.root_dir = os.path.dirname(os.path.abspath(__file__))
		req = open(kwargs["requirement"])
		for line in req:
			rows = line.strip().split()
			if len(rows) < 6:
				continue
			s_fname, fname, version, key, passwd, ftype = rows[:6]
			if ftype not in ['db', 'source']:
				continue
			self.files[ftype][s_fname] = (fname, key, passwd, version)
		req.close()

	def install(self):
		if not self.nodata:
			for db in self.files['db']:
				self.download(db, 'db')
		if not self.nosoft:
			if 'bwa' in self.files['source']:
				self.install_bwa(versions=self.files['source']['bwa'][-1])
			if 'SOAPnuke' in self.files['source']:
				self.install_sopenuke()
			if 'samtools' in self.files['source']:
				self.install_samtools(versions=self.files['source']['samtools'][-1])
			if 'picard' in self.files['source']:
				self.install_picard()
			if 'GATK3' in self.files['source']:
				self.install_gatk()
			if 'snpEff' in self.files['source']:
				self.install_snpeff()
		for rh_script in os.listdir(os.path.join(self.root_dir, "tools")):
			if os.path.exists(os.path.join(self.prefix, "tools", rh_script)):
				continue
			if os.path.isfile(rh_script):
				shutil.copy(rh_script, os.path.join(self.prefix, "tools"))
			elif os.path.isdir(rh_script):
				shutil.copytree(rh_script, os.path.join(self.prefix, "tools"), symlinks=True)
			else:
				continue

	def _download(self, outdir, filename, dbkey, dbword):
		out_fname = os.path.basename(outdir)
		fname = os.path.join(outdir, filename)
		logger.info("* downloading %s to %s " % (out_fname, outdir))
		url = "/".join([self._download_url, dbkey, filename, dbword])
		retries = 0
		while retries < self.max_retries:
			cmd = ["wget", "--continue", "--quiet", "-O", fname, url]
			retcode = call(cmd)
			if retcode == 0:
				break
			logger.info("wget failed with non-zero exit code %s. Retrying" % retcode)
			retries += 1
			if retries >= self.max_retries:
				raise ValueError("Failed to download with wget")
			time.sleep(10)
		if filename.endswith('tar.gz'):
			tar = tarfile.open(fname, "r:gz")
			file_names = tar.getnames()
			for file_name in file_names:
				tar.extract(file_name, outdir)
			tar.close()
			os.remove(fname)

	def download(self, name, ftype):
		if name not in self.files[ftype]:
			raise ValueError("not support %s yet !" % name)
		logger.info("* installing %s files : %s " % (ftype, name))
		outdir = os.path.join(self.prefix, ftype, name)
		if not os.path.isdir(outdir):
			os.makedirs(outdir)
		raw_name, key, passwd, version = self.files[ftype][name]
		self._download(outdir, raw_name, key, passwd)

	@staticmethod
	def is_tool(soft, version=None):
		try:
			s_i, s_o = Popen(soft.strip().split(), stdout=PIPE, stderr=PIPE).communicate()
			if s_o and version is not None:
				for i in s_o.strip().split("\n"):
					rows = i.strip().split()
					if not len(rows) > 1:
						continue
					if rows[0].lower().startswith("version") and LooseVersion(rows[1]) < LooseVersion(version):
						return False
			if s_i:
				return s_i
		except OSError as e:
			if e.errno == os.errno.ENOENT:
				return False
			else:
				raise OSError(e)
		return True

	def install_bwa(self, versions=None):
		newpath = os.path.join(self.prefix, "tools", "bwa")
		if self.is_tool("bwa", versions):
			old_path = self.is_tool("which bwa").strip()
			if old_path and os.path.isfile(old_path):
				mksymlink(old_path, newpath)
		else:
			self.download('bwa', ftype="source")
			outdir = os.path.join(self.prefix, 'source', 'bwa')
			os.chdir(outdir)
			_ = self.is_tool("make")
			os.chdir(self.root_dir)
			mksymlink(os.path.join(outdir, 'bwa'), newpath)

	def install_samtools(self, versions=None):
		newpath = os.path.join(self.prefix, "tools", "samtools")
		if self.is_tool("samtools", versions):
			old_path = self.is_tool("which samtools").strip()
			if old_path and os.path.isfile(old_path):
				mksymlink(old_path, newpath)
		else:
			self.download('samtools', ftype="source")
			outdir = os.path.join(self.prefix, 'source', 'samtools')
			os.chdir(outdir)
			_ = self.is_tool("./configure --prefix=%s " % outdir)
			_ = self.is_tool("make")
			_ = self.is_tool("make install")
			os.chdir(self.root_dir)
			mksymlink(os.path.join(outdir, 'bin', 'samtools'), newpath)

	def install_picard(self):
		newpath = os.path.join(self.prefix, "tools", "picard")
		self.download('picard', ftype="source")
		mksymlink(os.path.join(self.prefix, 'source', 'picard'), newpath)

	def install_snpeff(self):
		newpath = os.path.join(self.prefix, "tools", "snpEff")
		self.download('snpEff', ftype="source")
		mksymlink(os.path.join(self.prefix, 'source', 'snpEff'), newpath)

	def install_gatk(self):
		newpath = os.path.join(self.prefix, "tools", "GATK3")
		self.download('GATK3', ftype="source")
		mksymlink(os.path.join(self.prefix, 'source', 'GATK3'), newpath)

	def install_sopenuke(self):
		newpath = os.path.join(self.prefix, "tools", "SOAPnuke")
		self.download('SOAPnuke', ftype="source")
		mksymlink(os.path.join(self.prefix, 'source', 'SOAPnuke'), newpath)


def main():
	usage = 'Usage: %prog [-h] [--version] --samplelist [sample list] --config [config file] [options]'
	description = 'scripts for install RheaChip'
	author = 'joeyhwong@hotmail.com'
	version = '0.1.0'
	parser = OptionParser(usage=usage, version=version, description=description, epilog=author)
	parser.add_option('-i', '--source', metavar='FILE', dest='sourcelist', default=None)
	parser.add_option('--prefix', metavar='DIR', dest="prefix", help='install dir', default=Rhea_Chip_Home)
	parser.add_option('--nodata', help='Do not install data dependencies', dest='nodata', default=False,
	                  action="store_true")
	parser.add_option('--nosoft', help='Do not install software dependencies', dest='nosoft', default=False,
	                  action="store_true")
	options, args = parser.parse_args()
	workdir = os.path.dirname(os.path.abspath(__file__))
	sourcelist = options.sourcelist or os.path.join(workdir, "source.list")
	prefix = options.prefix or workdir
	if not os.path.isfile(sourcelist):
		parser.print_help()
		return 'Expected source lists lost !!!'
	jobs = Install(prefix=prefix, requirement=sourcelist, nodata=options.nodata, nosoft=options.nosoft)
	jobs.install()


if __name__ == '__main__':
	sys.exit(main())
