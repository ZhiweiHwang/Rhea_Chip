#! /usr/bin/env python
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！
import numpy as np


class ViterbiTraining(object):
	def __init__(self, data, ploid, best_probability, trials):
		self.data = data
		self.ploid = float(ploid)
		self.best_probability = float(best_probability)
		self.trials = int(trials)
		self.ini_P = {0: 0.01, 1: 0.03, 2: 0.9, 3: 0.03, 4: 0.02, 5: 0.01}
		self.end_P = {0: 0.01, 1: 0.03, 2: 0.9, 3: 0.03, 4: 0.02, 5: 0.01}
		self.min_trans_P = min(self.ini_P.values())
		self.CN = len(self.ini_P.keys())

	def mlEstimate(self, copy_est):
		trans_p = np.array([[self.min_trans_P, ] * self.CN] * self.CN)
		for i in xrange(len(copy_est) - 1):
			trans_p[copy_est[i]][copy_est[i + 1]] += 1
		sum_p = trans_p.sum(axis=1)
		for i in xrange(self.CN):
			if (sum_p[i] == self.CN * self.min_trans_P) or (sum_p[i] < 1.0):
				trans_p[i][i] += 1
				sum_p[i] += 1
			trans_p[i] /= sum_p[i]
		return trans_p

	def emi(self, copyn, dep):
		if copyn == 0:
			return np.exp(-0.35 * dep) / 3.3863 + 0.000001
		nbinom_p = self.best_probability
		ans = self.best_probability
		dd = int(self.ploid * dep / copyn + 0.5)
		if (copyn == self.CN - 1) and (dep > 30 * copyn):
			dd = 60
		for i in xrange(self.trials - 1):
			ans *= nbinom_p * (dd + i + 1.0) / (i + 1.0)
		for i in xrange(dd):
			ans *= (1 - nbinom_p)
		return ans + 0.000001

	def viterbi(self, trans_p):
		trans_p = np.log10(trans_p + 0.000001)
		vi = np.zeros([len(self.data), self.CN, 2])
		copy_est = np.zeros(len(self.data))
		for i in range(self.CN):
			vi[0][i][0] = np.log10(self.emi(i, self.data[0]) + 0.000001) + np.log10(self.ini_P[i])
			vi[0][i][0] = -1
		for i in xrange(1, len(self.data)):
			for j in range(self.CN):
				maxvt = - np.sys.maxint
				maxk = 0
				for k in range(self.CN):
					vt = vi[i - 1][k][0] + trans_p[k][j]
					if vt > maxvt:
						maxvt = vt
						maxk = k
				vi[i][j][0] = maxvt + np.log10(self.emi(j, self.data[i]) + 0.000001)
				vi[i][j][1] = maxk
		end_c = 0
		end_vt = - np.sys.maxint
		for i in range(self.CN):
			if vi[-1][i][0] > end_vt:
				end_vt = vi[-1][i][0]
				end_c = i
		copy_est[-1] = end_c
		for i in range(len(self.data) - 1, 0, -1):
			copy_est[i - 1] = int(vi[i][copy_est[i]][1])
		return copy_est

	def train(self, copy_est=None):
		if copy_est is None:
			mdep = np.mean(self.data) * 2.0 / self.ploid
			copy_est = list()
			for i in self.data:
				dr = i / mdep
				if dr < 0.1:
					copy_est.append(0)
				elif dr < 0.7:
					copy_est.append(1)
				elif dr <= 1.35:
					copy_est.append(2)
				elif dr <= 1.75:
					copy_est.append(3)
				elif dr <= 2.25:
					copy_est.append(4)
				else:
					copy_est.append(5)
		trans_p = self.mlEstimate(copy_est)
		indi = 1
		while indi:
			n_copy_est = self.viterbi(trans_p)
			trans_p = self.mlEstimate(n_copy_est)
			indi = np.count_nonzero(np.array(n_copy_est) != copy_est)
			copy_est = np.copy(n_copy_est)
		return copy_est

	def forward(self, tran):
		fs = np.zeros(len(self.data))
		fw = np.zeros((len(self.data), self.CN))
		for i in xrange(self.CN):
			fw[0][i] = self.ini_P[i] * self.emi(i, self.data[0])
			fs[0] += fw[0][i]
		fw[0] /= fs[0]
		for i in xrange(1, len(self.data)):
			for j in xrange(self.CN):
				fw[i][j] = sum([fw[i - 1][k] * tran[k][j] for k in xrange(self.CN)])
				fw[i][j] *= self.emi(j, self.data[i])
				fs[i] += fw[i][j]
			fw[i] /= fs[i]
		ans = np.log10(sum([fw[-1][j] * self.end_P[j] for j in xrange(self.CN)])) + sum(np.log10(fs))
		return ans, fs, fw

	def backward(self, tran):
		bw = np.zeros((len(self.data), self.CN))
		bs = np.zeros(len(self.data))
		for i in xrange(self.CN):
			bw[-1, i] = self.end_P[i]
		bs[-1] = sum(bw[-1, :])
		bw[-1] /= bs[-1]
		for i in xrange(len(self.data) - 2, -1, -1):
			for j in xrange(self.CN):
				bw[i][j] = sum([tran[j][k] * self.emi(k, self.data[i + 1]) * bw[i + 1][k] for k in xrange(self.CN)])
				bs[i] += bw[i][j]
			bw[i] /= bs[i]
		ans = np.log10(sum([bw[0][j] * self.ini_P[j] * self.emi(j, self.data[0]) for j in xrange(self.CN)]))
		ans += sum(np.log10(bs))
		return ans, bs, bw

	def posterior_decoding(self, trans, plody=2.0):
		fwp, fs, fw = self.forward(trans)
		bwp, bs, bw = self.backward(trans)
		avp = (fwp + bwp) / 2.0
		maxcn = float(plody)
		cne = list()
		for i in xrange(len(self.data)):
			logpp = sum(np.log10(fs[:i + 1])) + sum(np.log10(bs[i:]))
			maxp = 0
			for j in xrange(self.CN):
				logpp_j = pow(10, logpp + np.log10(fw[i, j]) + np.log10(bw[i, j]) - avp)
				if logpp_j > maxp:
					maxp = logpp_j
					maxcn = j
			cne.append([maxp, maxcn])
		return cne
