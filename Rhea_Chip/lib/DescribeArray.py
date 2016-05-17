#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！
import numpy as np


class DescribeArray(object):
	def __init__(self, numberlist, col=None, dtype=object):
		self.array = np.array(numberlist, dtype=dtype)
		if col is not None:
			self.average = self.array[:, int(col)].mean()
			self.median = np.median(self.array[:, int(col)])
			self.max = self.array[:, int(col)].max()
			self.min = self.array[:, int(col)].min()
		else:
			self.average = np.mean(self.array)
			self.median = np.median(self.array)
			self.max = np.max(self.array)
			self.min = np.min(self.array)

	def __str__(self):
		return self.average, self.median, self.max, self.min

	def get_frequece(self, thresh, col=None):
		if col is not None:
			radio = round(np.count_nonzero(self.array[:, int(col)] >= float(thresh)) / float(len(self.array)), 4) * 100
			return "%4.2f%%" % radio
		else:
			radio = round(np.count_nonzero(self.array >= float(thresh)) / float(len(self.array)), 4) * 100
			return "%4.2f%%" % radio
