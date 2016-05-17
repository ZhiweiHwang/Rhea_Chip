#! /usr/bin/env python
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！
import numpy as np


class NegativeBinomial(object):
	def __init__(self, data, runtime=290, train=4):
		data = [int(float(i) + 0.5) for i in data]
		self.train = int(train) or 4
		self.runtimes = int(runtime) or 290
		self.max = max(data)
		self.mindevi = self.max + 0.1
		self.data = [float(i) / len(data) for i in np.bincount(data)]
		self.devi = np.zeros(11)
		self.probabilities = np.arange(0, 1.1, 0.1)
		self.best_probability = 0.0
		self.trials = 0

	def nbinom_fit(self):
		bp = 0
		for runs in range(2, self.runtimes):
			for _ in range(self.train):
				bp = self.nbfitp(runs)
			if self.devi[bp] < self.mindevi:
				self.mindevi = round(self.devi[bp], 6)
				self.trials = runs
				self.best_probability = self.probabilities[bp]

	def nbfitp(self, r):
		p = np.arange(0, 1.1, 0.1)
		devi = np.zeros(11)
		bp = 0
		gap = p[1] - p[0]
		mindevi = self.max + 0.1
		for i in xrange(11):
			for j in xrange(self.max + 1):
				deviation = self.data[j] - self.nbinop(j, r, p[i])
				devi[i] += deviation * deviation
			if devi[i] < mindevi:
				mindevi = devi[i]
				bp = i
		if bp == 0:
			ans = 0
		elif bp == 10:
			p[0] = p[9]
			ans = 10
		elif devi[bp - 1] < devi[bp + 1]:
			ans = 10
			p[0] = p[bp - 1]
		else:
			ans = 0
			p[0] = p[bp]
		gap /= 10.0
		self.probabilities = [p[0] + gap * i for i in range(11)]
		devi[ans] = mindevi
		self.devi = devi
		return ans

	@staticmethod
	def nbinop(k, r, p):
		ans = p
		for i in range(int(r) - 1):
			ans *= p * (int(k) + i + 1.0) / (i + 1.0)
		ans *= pow(1.0 - p, int(k))
		return ans
