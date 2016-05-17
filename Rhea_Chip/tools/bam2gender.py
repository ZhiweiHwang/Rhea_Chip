#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
get gender from bamfile
"""

from __future__ import absolute_import, unicode_literals, division
import sys
import os
import pysam
import re
from toolz import partial, pipe
from toolz.curried import map
import numpy as np


def average(sequence):
	if not len(sequence):
		return 0.0
	try:
		return sequence.mean()
	except AttributeError:
		return sum(sequence) / len(sequence)


class Predict(object):
	"""
	gender prediction based on the depth
	"""

	def __init__(self, bamfile, outdir):
		self.bamfile = bamfile
		stat = self.indexbamfile()
		self.outdir = outdir
		assert self.bamfile and self.outdir and stat, "Input error"
		self._bam = pysam.Samfile(bamfile)
		self._prealloc_func = partial(np.zeros, dtype=np.int)
		self.fake_bed_rows = [("chrX", 1, 59373566), ("chrY", 69362, 11375310)]
		self.sequence = pipe(self.fake_bed_rows,
		                     map(lambda interval: self.depthreader(*interval)),
		                     map(average)
		                     )
		self.x_coverage, self.y_coverage = list(self.sequence)
		self.sex = self.predict_gender()

	def indexbamfile(self):
		bamdir, bamfile = os.path.split(self.bamfile)
		if os.path.isfile(self.bamfile + '.bai') or os.path.isfile(re.sub("bam$", 'bai', self.bamfile)):
			return True
		if not os.access(bamdir, os.W_OK):
			try:
				os.symlink(self.bamfile, os.path.join(self.outdir, bamfile))
				self.bamfile = os.path.join(self.outdir, bamfile)
			except Exception:
				return False
		try:
			_ = pysam.index(self.bamfile, force=True)
			return True
		except Exception:
			return False

	def depthreader(self, contig, start, end):
		pysam_start = start - 1
		pysam_contig = str(contig)
		if pysam_start < 0:
			raise ValueError("Start position must be > 0, not %d" % start)
		read_depths = self._prealloc_func(end - pysam_start)
		for col in self._bam.pileup(pysam_contig, pysam_start, end, truncate=True):
			read_depths[col.pos - pysam_start] = col.n
		return read_depths

	def predict_gender(self):
		if self.x_coverage == 0:
			return 'unknown'
		elif (self.y_coverage > 0) and (self.x_coverage / self.y_coverage < 10):
			return 'male'
		else:
			return 'female'


def main():
	usage = 'python %s <bamfile> <outfile>' % sys.argv[0]
	if len(sys.argv) != 3:
		return usage
	bamfile = os.path.realpath(sys.argv[1])
	outfile = os.path.realpath(sys.argv[2])
	outdir, name = os.path.split(outfile)
	out = open(outfile, 'w')
	predict = Predict(bamfile, outdir)
	out.writelines("\t".join([name, str(predict.sex), str(predict.x_coverage), str(predict.y_coverage)]) + '\n')
	out.close()
	return 0


if __name__ == '__main__':
	sys.exit(main())
