#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！
import numpy
import matplotlib
matplotlib.use('Agg')
import multiprocessing
import os
import pysam
import re
from collections import defaultdict
from Rhea_Chip.lib.hgvs import HGVS
from Rhea_Chip.lib.DescribeArray import DescribeArray
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def cpu_count():
	return multiprocessing.cpu_count()


def _exon_vaule(exons):
	v = re.compile("[0-9]+$").match(exons)
	if v:
		return int(v.group())
	else:
		return sum([ord(i) for i in exons])


def graph(depthfile, messages, outdir=os.getcwd(), prefix="Test", uncoverthreshold=5):
	depth = pysam.TabixFile(depthfile)
	plot_data = defaultdict(list)
	genes = str
	trans = str
	for exon_message in sorted(messages, key=lambda x: _exon_vaule(x[2])):
		if len(exon_message) != 6:
			continue
		trans, genes, exon_name, refchrom, refstart, refend = exon_message
		try:
			chrom = str(refchrom)
			if chrom.startswith("chrM"):
				chrom = 'chrM_NC_012920.1'
			start = int(refstart)
			end = int(refend)
		except Exception:
			continue
		if chrom not in depth.contigs:
			continue
		array = list()
		for depth_m in depth.fetch(chrom, start - 1, end):
			rows = depth_m.strip().split('\t')
			if len(rows) < 3:
				raise IOError("depth file format Error !")
			d = rows[-1]
			array.append(int(d))
		if not len(array):
			array = [0] * (end - start + 1)
		depth_message = DescribeArray(array)
		plot_data["averagedepth"].append(depth_message.average)
		plot_data["mediandepth"].append(depth_message.median)
		coverate = float(depth_message.get_frequece(uncoverthreshold).strip('%'))
		plot_data["coverate"].append(coverate)
		plot_data["x_axis"].append(exon_name)
	depth.close()
	plot(plot_data, outdir, prefix, genes, trans)


def plot(plot_data, outdir, prefix, genes, trans):
	prefixname = ".".join([prefix, genes, trans])
	output = os.path.join(outdir, prefixname + '.png')
	size = max(len(plot_data["x_axis"]), 8)
	fig, axes = plt.subplots(2, 1)
	plt.suptitle(prefixname, fontsize=18)
	depth_ax, cover_ax = axes
	depth_ax.get_xaxis().set_major_locator(MaxNLocator(nbins=len(plot_data["x_axis"]), prune='upper'))
	artists1 = depth_ax.bar(numpy.arange(len(plot_data["x_axis"])) + 0.2, numpy.array(plot_data["averagedepth"]),
	                        width=0.4, label="Mean depth", align="center", color='c')
	artists2 = depth_ax.bar(numpy.arange(len(plot_data["x_axis"])) + 0.6, numpy.array(plot_data["mediandepth"]),
	                        width=0.4, color="pink", label="Median depth", align="center")
	lefts = [item.get_x() for item in artists1]
	rights = [item.get_x() + item.get_width() for item in artists2]
	depth_ax.set_xlim([min(lefts), max(rights)])
	depth_ax.set_title('Depths on all exon regions of %s' % trans, fontsize=15)
	depth_ax.set_ylabel("Depth (X)", fontsize=13)
	labels = numpy.array(plot_data["x_axis"])
	depth_ax.set_xticklabels(labels, rotation=45, rotation_mode="anchor", size=min(180.0 / size, 13))
	depth_ax.legend(loc='best', prop={'size': 10})
	cover_ax.get_xaxis().set_major_locator(MaxNLocator(nbins=len(plot_data["x_axis"]), prune='upper'))
	cover_ax.plot(numpy.arange(len(plot_data["x_axis"])) + 0.6, numpy.array(plot_data["coverate"]),
	              'yo-', label="Coverage radio")
	cover_ax.set_title('Coverages on all exon regions of %s' % trans, fontsize=15, verticalalignment='bottom')
	cover_ax.set_ylabel("Radio (%)", fontsize=13)
	cover_ax.set_ylim(0, 105)
	cover_ax.set_xticklabels(labels, rotation=45, rotation_mode="anchor", size=min(180.0 / size, 13))
	cover_ax.get_xaxis().set_visible(False)
	cover_ax.legend(loc='best', prop={'size': 10})
	plt.savefig(output, bbox_inches='tight', dpi=240, figsize=(size * 2, size))
	plt.close()


def exon_graph(depth, refseq, prefix, outdir=os.getcwd(), trans=None, genes=None, threads=0):
	beds = HGVS(refseq, trans=trans, genes=genes)
	messages = beds.get_exons()
	threads = cpu_count() if not int(threads) else int(threads)
	pool = multiprocessing.Pool(processes=threads)
	for transcript, trans_messages in messages.iteritems():
		pool.apply_async(graph, (depth, trans_messages, outdir, prefix))
	pool.close()
	pool.join()
