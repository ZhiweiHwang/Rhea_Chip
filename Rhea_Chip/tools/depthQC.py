#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！
import os
import sys
import re
import pysam
import Rhea_Chip
from collections import defaultdict
from glob import glob
from Rhea_Chip.lib.CreateFolder import FolderMaker
from Rhea_Chip.lib.CreatIndex import CreatIndex
from Rhea_Chip.lib.bedAnno import BedAnno
from Rhea_Chip.lib.DescribeArray import DescribeArray
from Rhea_Chip.lib.ExonGraph import exon_graph
from Rhea_Chip.tools import _chrom_valued, formact_number_list_to_range
from smart_open import smart_open
from optparse import OptionParser, OptionGroup


class BamDepth(object):
	def __init__(self, bamfile, reference):
		bam = CreatIndex(bamfile)
		try:
			self._reader = pysam.AlignmentFile(bam.files, 'rb')
			self.name = str(self._reader.header['RG'][0]['SM'])
			self.reference = pysam.FastaFile(reference)
		except ValueError:
			raise ValueError('File header is empty, input bam / reference error !')

	def __del__(self):
		self._reader.close()
		self.reference.close()

	@staticmethod
	def count_gc(bases):
		try:
			gc_base = len(re.findall("[GC]", bases))
			total_base = len(re.findall("[GCTA]", bases))
			return gc_base * 100.0 / total_base
		except ZeroDivisionError:
			return 0.0

	def depths(self, bed, prefix=None, read_filter='all', qual_threshold=15, count_threshold=None, uncover_threshold=5):
		threshold_type = type(count_threshold)
		if threshold_type is list or threshold_type is tuple:
			count_threshold = set(count_threshold)
		elif count_threshold is set:
			pass
		else:
			count_threshold = {1, 4, 10, 20, 30, 100}
		uncover_threshold = max(int(uncover_threshold), 1)
		count_threshold.add(uncover_threshold)
		count_threshold = sorted(count_threshold, key=int)
		outdir, name = os.path.split(os.path.abspath(prefix)) if prefix else (os.getcwd(), 'Rhea_Chip')
		depth = os.path.join(outdir, '%s.depth.tsv' % name)
		bedstat = os.path.join(outdir, '%s.bed.stat' % name)
		stats = os.path.join(outdir, '%s.stat' % name)
		chromstat = os.path.join(outdir, '%s.chrom.stat' % name)
		uncover = os.path.join(outdir, '%s.uncover.bed' % name)
		_depth = smart_open(depth, 'w')
		_bedstat = smart_open(bedstat, 'w')
		_stats = smart_open(stats, 'w')
		_chromstat = smart_open(chromstat, 'w')
		# _depth.write("#Chrom\tPos\tRef\tCov_A\tCov_C\tCov_G\tCov_T\tDepth\tWinGC_%s\n" % str(gcwin))
		_depth.write("#Chrom\tPos\tRef\tCov_A\tCov_C\tCov_G\tCov_T\tDepth\n")
		_bedstat.write("#Chr\tStart\tStop\tAverage\tMedian\tMax\tMin\n")
		_stats.write("##A Simple introduction about %s\n" % self.name)
		_chromstat.write("#Chr\tAverage\tMedian\tMax\tMin\t")
		_chromstat.write("\t".join(["Coverage (>=%sX)" % str(key) for key in count_threshold]) + '\n')
		chroms = defaultdict(dict)
		rangedict = defaultdict(int)
		regiondict = defaultdict(int)
		total_base = 0
		region_num = 0
		with smart_open(bed) as regions:
			for region in regions:
				rows = region.strip().split()
				if len(rows) < 3:
					continue
				try:
					chrom = rows[0]
					start = max(int(rows[1]) - 1, 0)
					stop = min(int(rows[2]) + 1, self.reference.get_reference_length(chrom))
				except Exception:
					continue
				cov_a, cov_c, cov_g, cov_t = self._reader.count_coverage(chrom, start, stop, read_callback=read_filter,
				                                                         quality_threshold=int(qual_threshold))
				bases = self.reference.fetch(chrom, start, stop).upper()
				reg = list()
				chrom = "chr" + re.sub("^chr", "", chrom)
				if chrom.startswith("chrM"):
					chrom = 'chrM_NC_012920.1'
				for n in xrange(start, stop):
					offset = n - start
					dep = [bases[offset], cov_a[offset], cov_c[offset], cov_g[offset], cov_t[offset]]
					base_depth = sum(dep[1:])
					dep.append(base_depth)
					# gc_radio = round(self.count_gc(bases[offset:offset + 201]), 4)
					# dep.append(gc_radio)
					chroms[chrom][n] = dep
					for num in count_threshold:
						num = int(num)
						if base_depth >= num:
							rangedict[num] += 1
					total_base += 1
					reg.append(base_depth)
				region_num += 1
				array = DescribeArray(reg)
				averages = str(round(array.average, 2))
				mediandepth = str(round(array.median, 2))
				maxdepth = str(round(array.max, 2))
				mindepth = str(round(array.min, 2))
				_bedstat.write("\t".join([chrom, rows[1], rows[2], averages, mediandepth, maxdepth, mindepth]) + '\n')
				for num in count_threshold:
					num = int(num)
					if array.average >= num:
						regiondict[num] += 1
		uncover_range = list()
		n_array = list()
		for chrom in sorted(chroms.keys(), key=lambda x: _chrom_valued(x)):
			dep = sorted(chroms[chrom].iteritems(), key=lambda x: int(x[0]))
			array = list()
			for p, d in dep:
				_depth.write("\t".join([chrom, str(p), '\t'.join([str(i) for i in d])]) + '\n')
				array.append(d)
				n_array.append(d[5])
			array = DescribeArray(array, col=5)
			averages = str(round(array.average, 2))
			mediandepth = str(round(array.median, 2))
			maxdepth = str(round(array.max, 2))
			mindepth = str(round(array.min, 2))
			chromcover = [array.get_frequece(thre, col=5) for thre in count_threshold]
			_chromstat.write("\t".join([chrom, averages, mediandepth, maxdepth, mindepth] + chromcover) + '\n')
			if read_filter == 'all':
				uncover_bases = [int(p) for p, d in dep if d[5] < uncover_threshold]
				uncover_range.extend(formact_number_list_to_range(uncover_bases, tag=chrom))
		if read_filter == 'all':
			uncoverout = smart_open(uncover, 'w')
			uncoverout.write("#Chr\tStart\tStop\n")
			uncoverout.writelines(uncover_range)
			uncoverout.close()
		array = DescribeArray(n_array)
		_stats.write("Average depth : %.2f\n" % array.average)
		_stats.write("Median depth : %.2f\n" % array.median)
		_stats.write("Max depth : %.2f\n" % array.max)
		_stats.write("Min depth : %.2f\n" % array.min)
		for number in count_threshold:
			number = int(number)
			_stats.write("Coverage (>={0:d}X) : {1:.2f}%\n".format(
				number, (float(rangedict[number]) / total_base) * 100, 2))
		for number in count_threshold:
			number = int(number)
			_stats.write("Region Coverage (>={0:d}X) : {1:.2f}%\n".format(
				number, (float(regiondict[number]) / region_num) * 100, 2))
		_depth.close()
		_bedstat.close()
		_stats.close()
		_chromstat.close()
		dep_f = CreatIndex(depth)
		_ = dep_f.check_index(seq_col=0, start_col=1, end_col=1)


def main():
	usage = 'Usage: %prog [-h] [--version] --target [target bed] -b [bam file] -r [ref seq] [options]'
	description = Rhea_Chip.__description__
	author = Rhea_Chip.__author__
	version = Rhea_Chip.__version__
	parser = OptionParser(usage=usage, version=version, description=description, epilog=author)
	expect = OptionGroup(parser, 'Expected arguments', 'Caution: These parameters are necessary for depthQC.')
	expect.add_option('-b', '--bam', metavar='FILE', dest='bam', help='bam file')
	expect.add_option('-r', '--reference', metavar='FILE', dest='reference', help='humman reference')
	expect.add_option('--target', metavar='FILE', dest='target', help='target bed region file')
	parser.add_option_group(expect)
	optinal = OptionGroup(parser, 'Optional arguments', 'Caution: If you do not set these parameters in addition,'
	                                                    ' depthQC will select the default value.')
	optinal.add_option('-o', '--outdir', dest='outdir', default=os.getcwd(),
	                   help='The output dir [ default: %s ]' % os.getcwd())
	optinal.add_option('-t', '--threads', dest='threads', default=0,
	                   help='The max threads [ default: <cpus> ], only use in plot exon graph')
	optinal.add_option('--flank', dest='flank', default=None, help='Flank bed region file [ default: None ]')
	optinal.add_option('--trans', dest='trans', default=None, help='Trans list file [ default: None ]')
	optinal.add_option('--genes', dest='genes', default=None, help='Genes list file [ default: None ]')
	optinal.add_option('-p', '--plot', help='Plot exon depth graph', dest='plot', default=False, action="store_true")
	parser.add_option_group(optinal)
	(options, args) = parser.parse_args()
	if not options.bam or not options.reference or not options.target:
		parser.print_help()
		return 'Expected arguments lost !!!'
	bam = os.path.abspath(options.bam)
	reference = os.path.abspath(options.reference)
	target = os.path.abspath(options.target)
	outdir = os.path.abspath(options.outdir)
	threads = int(options.threads)
	flank = os.path.abspath(options.flank) if options.flank else None
	trans = os.path.abspath(options.trans) if options.trans else None
	genes = os.path.abspath(options.genes) if options.genes else None
	plot = options.plot
	assert trans or genes, "Genes and Trans must set one ~"

	prepair_dir = FolderMaker(outdir)
	subdir = ["depthAnno", "depthAnno/Target"]
	if flank:
		subdir.append("depthAnno/Flank")
	if plot:
		prepair_dir.create_subdirectory(["graph_exon"])
	prepair_dir.create_subdirectory(subdir)
	depths = BamDepth(bam, reference)
	sample = depths.name
	depths.depths(target, prefix=os.path.join(outdir, "depthAnno/Target", sample + '.Rmdup'))
	depths.depths(target, prefix=os.path.join(outdir, "depthAnno/Target", sample + '.Raw'), read_filter="nofilter")
	if flank:
		depths.depths(flank, prefix=os.path.join(outdir, "depthAnno/Flank", sample + '.Rmdup'))
		depths.depths(target, prefix=os.path.join(outdir, "depthAnno/Flank", sample + '.Raw'), read_filter="nofilter")
	depths.__del__()
	uncovers = glob(os.path.join(outdir, 'depthAnno', "*", sample + ".Rmdup.uncover.bed"))
	bedanno = BedAnno(reference=reference, trans=trans, genes=genes)
	for bed in uncovers:
		file_out = re.sub("\.bed$", ".anno", bed)
		bedanno.bedanno(bed, fileout=file_out)
	if plot:
		depth_f = os.path.join(outdir, "depthAnno/Target", sample + '.Rmdup.depth.tsv.gz')
		exon_graph(depth=depth_f, refseq=bedanno.refdb, prefix=sample,
		           outdir=os.path.join(outdir, "graph_exon"), trans=trans,
		           genes=genes, threads=threads)
	bedanno.__del__()


if __name__ == '__main__':
	sys.exit(main())
