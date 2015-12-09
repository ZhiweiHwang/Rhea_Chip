#! /usr/bin/env python
# -*- coding: utf-8 -*-
# File Name: depthQC.py

"""
QC for bam file
"""

from collections import defaultdict
from optparse import OptionParser
import glob
import itertools
import json
import matplotlib
import numpy
import os
import pysam
import re
import sys
import warnings

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def custom_formatwarning(msg, *a):
	# ignore everything except the message
	return "[Warning]" + str(msg) + '\n'


warnings.formatwarning = custom_formatwarning


class ConfigPhrase(object):
	def __init__(self, filein, fileout=None):
		"""
		Phrase config file
		"""
		self.Rheadir = os.getcwd()
		self.raw_config_file = filein
		self.config_save_file = fileout
		self.config = self.config_load
		self.tools = dict()
		self.config_stat = self.check_config

	@property
	def config_load(self):
		"""
		loadding all the personalized settings from json file and return a dict
		and store all personalized settings to a file in json format
		"""
		with open(self.raw_config_file) as config_file:
			config = json.load(config_file)
		if self.config_save_file is not None:
			with open(self.config_save_file, 'w') as f:
				json.dump(config, f, indent=4)
		try:
			self.Rheadir = config['Rhea_Chip_dir']['PATH']
		except ValueError:
			pass
		return config

	def get_real_path(self, indexpath):
		realpath = os.path.realpath(os.path.join(self.Rheadir, indexpath))
		return os.path.exists(realpath), realpath

	@property
	def check_config(self):
		"""
		Check that all required parameters are complete
		"""
		root_check = os.path.exists(self.Rheadir)
		try:
			for otherdata in ('dbsnp', 'Refseq'):
				temp_check, realpath = self.get_real_path(self.config['Datasets'][otherdata]['PATH'])
				root_check *= temp_check
				self.tools[otherdata] = realpath
		except ValueError:
			root_check = 0

		self.tools["Trans"] = set()
		self.tools["Gene"] = set()
		try:
			for k in ("Trans", "Gene"):
				temp_check, realpath = self.get_real_path(self.config["Personal Analysing"]["%s List" % k]['PATH'])
				root_check *= temp_check
				if temp_check:
					with open(realpath) as lists:
						for n in lists:
							if k is "Trans":
								n = re.sub('\.[^\.]+$', '', n)
							self.tools[k].add(n.strip())
		except ValueError:
			pass

		try:
			bedfile = self.config["Personal Analysing"]["Bed"]["Target bed"]
			temp_check_target, self.tools["Target bed"] = self.get_real_path(bedfile)
		except ValueError:
			pass

		try:
			bedfile = self.config["Personal Analysing"]["Bed"]["Flank bed"]
			temp_check_flank, self.tools["Flank bed"] = self.get_real_path(bedfile)
		except ValueError:
			pass

		try:
			if not self.tools["Flank bed"]:
				if not self.tools["Target bed"]:
					self.tools["Target bed"] = self.bed_region_take()
				self.tools["Flank bed"] = self.extend_bed_region(os.path.join(temppath, "flank.bed"))
			elif not self.tools["Target bed"]:
				self.tools["Target bed"] = self.tools["Flank bed"]
		except Exception:
			root_check = 0

		return root_check

	def bed_region_take(self):
		target_bed = open(os.path.join(temppath, "target.bed"), 'w')
		with open(self.tools['Refseq']) as refbed:
			for lines in refbed:
				rows = lines.strip().split('\t')
				if len(rows) < 9:
					continue
				trans = re.sub('\.[^\.]+$', '', rows[0])
				genes = rows[2]
				chrom = rows[3]
				start = rows[4]
				end = rows[5]
				region = rows[8]
				if region.startswith('IVS') or not start.isdigit() or not end.isdigit() or not chrom.startswith("chr"):
					continue
				if len(self.tools["Trans"]) and trans not in self.tools["Trans"]:
					continue
				if len(self.tools["Gene"]) and genes not in self.tools["Gene"]:
					continue
				if not len(self.tools["Trans"]):
					self.tools["Trans"].add(trans)
				if not len(self.tools["Gene"]):
					self.tools["Gene"].add(genes)
				target_bed.write("\t".join([chrom, start, end]) + '\n')
		target_bed.close()
		return target_bed

	def extend_bed_region(self, fileout=None, flank=100):
		from itertools import chain
		flatten = chain.from_iterable
		left, right, check, x = 1, -1, 0, 0
		target_bed = defaultdict(list)
		fileouts = open(fileout, 'w') if fileout is not None else sys.stdout
		with open(self.tools["Target bed"]) as regions:
			for region in regions:
				rows = region.strip().split()
				if len(rows) < 3:
					continue
				try:
					chrom = str(regions[0])
					start = int(regions[1])
					end = int(regions[2])
					target_bed[chrom].append([max(start - flank, 0), max(end + flank, 0)])
				except:
					fileouts.write(region)
					continue
		for chrom in sort_bed_name_list(target_bed.keys()):
			for value, label in sorted(flatten(((start, left), (stop, right)) for start, stop in target_bed[chrom])):
				if check == 0:
					x = value
				check += label
				if check == 0:
					fileouts.write("\t".join([chrom, str(x), str(value)]) + '\n')
		fileouts.close()
		if fileout is not None:
			return fileout


class CreateFolder(object):
	"""
	a simple class use for makedir
	"""

	def __init__(self, rootdir):
		self.rootdir = self.folder_confirm(os.path.realpath(rootdir))

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


class BedAnno(object):
	def __init__(self, tools=None, configfile=None):
		if type(tools) != dict:
			if configfile is None:
				raise TypeError("BedAnno at least 1 positional argument (0 given)")
			else:
				baseconfig = ConfigPhrase(configfile)
				tools = baseconfig.tools
		self._refseq = pysam.TabixFile(tools['Refseq'])
		self.trans = tools["Trans"]
		self.genes = tools["Gene"]
		try:
			self._dbsnp = pysam.TabixFile(tools['dbsnp'])
		except Exception:
			self._dbsnp = None
		self.trans_region, self.target_anno = self.bedanno(tools["Target bed"])
		_, self.flank_anno = self.bedanno(tools["Flank bed"])

	def __del__(self):
		self._refseq.close()
		self._dbsnp.close()

	def bedanno(self, bedfile):
		transdict = defaultdict(list)
		annoresult = list()
		annoresult.append("#Transcript\tCodingRegion\tGenename\tChrom\tStart\tEnd\t"
		                  "Strand\tExon\tc.hgvs(start)\tc.hgvs(end)\n")
		with open(bedfile) as regions:
			for region in regions:
				rows = region.strip().split()
				if len(rows) < 3:
					continue
				try:
					chrom = str(rows[0])
					start = int(rows[1])
					end = int(rows[2])
				except Exception as err:
					warnings.warn("In line %s: %s" % (region.strip(), err))
					continue
				if chrom not in self._refseq.contigs:
					warnings.warn('Chromosome\'s name : %s is not in Refseq\'s contigs !!' % chrom)
					continue
				try:
					refannos = self._refseq.fetch(chrom, start, end)
				except Exception as err:
					warnings.warn("In line %s: %s" % (region.strip(), err))
					continue
				for refanno in refannos:
					refseq = refanno.strip().split('\t')
					if len(refseq) != 11:
						continue
					trans = refseq[0]
					codingregion = refseq[1]
					genes = refseq[2]
					refchrom = refseq[3]
					refstart = refseq[4]
					refend = refseq[5]
					strand = refseq[7]
					exon_name = refseq[8]
					cdsstart = refseq[9]
					cdsend = refseq[10]
					if len(self.trans) and re.sub('\.[^\.]+$', '', trans) not in self.trans:
						continue
					if len(self.genes) and genes not in self.genes:
						continue
					transdict[trans].append((genes, strand, exon_name, codingregion, refchrom, refstart, refend))
					start_hgvs = 'c.' + cdsstart
					end_hgvs = 'c.' + cdsend
					start_pos = refstart
					end_pos = refend
					if int(refend) > start > int(refstart):
						if (start - int(refstart) + 1) < (int(refend) - start + 1):
							start_hgvs = self.format_hgvs_name(cdsstart, start - int(refstart), exon_name, strand)
						else:
							_strand = '-' if strand == '+' else '+'
							start_hgvs = self.format_hgvs_name(cdsend, int(refend) - start, exon_name, _strand)
						start_pos = str(start)
					if int(refend) > end > int(refstart):
						if (end - int(refstart) + 1) < (int(refend) - end + 1):
							end_hgvs = self.format_hgvs_name(cdsstart, end - int(refstart), exon_name, strand)
						else:
							_strand = '-' if strand == '+' else '+'
							end_hgvs = self.format_hgvs_name(cdsend, int(refend) - end, exon_name, _strand)
						end_pos = str(end)
					annoresult.append("\t".join([trans, codingregion, genes, refchrom, start_pos, end_pos, strand,
					                             exon_name, start_hgvs, end_hgvs]) + '\n')
		return transdict, annoresult

	@staticmethod
	def format_hgvs_name(cds_pos, distance, exon, strand):
		if not len(cds_pos):
			return '-'
		if exon.startswith('IVS'):
			return 'c.' + cds_pos + strand + str(distance + 1)
		else:
			tag = ""
			if cds_pos.startswith('*'):
				tag, cds_pos = re.compile("(\*)(\d+)").match(cds_pos).groups()
			hgvs = int(cds_pos) + distance if strand == '+' else int(cds_pos) - distance
			return 'c.' + tag + str(hgvs)

	def uncoveranno(self, uncoverbed, fileout=None):
		f_out = open(fileout, 'w') if fileout is not None else sys.stdout
		_, tmp_anno_result = self.bedanno(uncoverbed)
		if self._dbsnp is None:
			f_out.writelines(tmp_anno_result)
		else:
			dbsnpanno = self._dbsnp
			for line in tmp_anno_result:
				line = line.strip()
				if line.startswith('#'):
					f_out.write(line + "\tMuts_in_dbsnp\n")
					continue
				rows = line.strip().split('\t')
				try:
					chrom = str(rows[3])
					start = int(rows[4])
					end = int(rows[5])
				except Exception as err:
					warnings.warn("In line %s: %s" % (line, err))
					f_out.write(line + "\t.\n")
					continue
				rawmut = set()
				for dbsnp_mut in dbsnpanno.fetch(chrom, start, end):
					mutchr, mutpos, mutid, mutref, mutcall = dbsnp_mut.split('\t')[:5]
					if mutid.startswith('rs'):
						rawmut.add(mutid)
					else:
						rawmut.add(":".join([mutchr, mutpos, mutref, mutcall]))
				if not len(rawmut):
					rawmut.add(".")
				f_out.write(line + "\t" + ";".join(rawmut) + "\n")
		f_out.close()
		if fileout is not None and os.path.isfile(fileout):
			os.remove(uncoverbed)
			return fileout
		else:
			return True


class BamDepth(object):
	def __init__(self, bamfile):
		try:
			self._bamreader = pysam.AlignmentFile(bamfile, 'rb')
			self.sample_name = str(self._bamreader.header['RG'][0]['SM'])
		except Exception as err:
			raise Exception(err)

	def __del__(self):
		self._bamreader.close()

	def countdepthcover(self, bed, fileout, stepper='all', count_threshold=None, uncover_threshold=4):
		tempcheck = type(count_threshold)
		if tempcheck is not set:
			if tempcheck is list or tempcheck is tuple:
				count_threshold = {int(i) for i in count_threshold}
			elif tempcheck is dict:
				count_threshold = {int(i) for i in count_threshold.keys()}
			else:
				count_threshold = {1, 4, 10, 20, 30, 100}
		uncover_threshold = max(int(uncover_threshold), 1)
		count_threshold.add(uncover_threshold)
		count_threshold = sorted(list(count_threshold))
		depths_out = open(fileout, 'w')
		coverout = open(fileout + '.stat', 'w')
		beddepth = open(fileout + '.depth', 'w')
		chromdepth = open(fileout + '.chrom.stat', 'w')
		depths_out.write("#Chr\tPos\tDepth\n")
		beddepth.write("#Chr\tStart\tEnd\tAverage\tMedian\tMax\tMin\n")
		chromdepth.write("#Chr\tAverage\tMedian\tMax\tMin")
		for key in count_threshold:
			chromdepth.write("\tCoverage (>=%sX)" % str(key))
		chromdepth.write('\n')
		chroms = defaultdict(dict)
		rangedict = defaultdict(int)
		regiondict = defaultdict(int)
		total_base = 0
		region_num = 0
		with open(bed) as regions:
			for region in regions:
				rows = region.strip().split('\t')
				depthlist = list()
				depths = defaultdict(int)
				if len(rows) < 3:
					continue
				try:
					chrom = str(rows[0])
					start = int(rows[1])
					end = int(rows[2])
				except Exception as err:
					warnings.warn("In line %s: %s" % (region.strip(), err))
					continue
				if chrom not in self._bamreader.references:
					warnings.warn('Chromosome\'s name : %s is not in bamfile\'s contigs !!' % chrom)
					continue
				try:
					pileups = self._bamreader.pileup(chrom, start - 1, end, stepper=stepper, truncate=True)
				except Exception as err:
					warnings.warn("In line %s: %s" % (region.strip(), err))
					continue
				for pileup in pileups:
					depths[pileup.pos + 1] = pileup.nsegments
				for pos in xrange(start, end + 1):
					depths_out.write("\t".join([chrom, str(pos), str(depths[pos])]) + '\n')
					chroms[chrom][pos] = depths[pos]
					depthlist.append(depths[pos])
					for num in count_threshold:
						num = int(num)
						if depths[pos] >= num:
							rangedict[num] += 1
				total_base += end - start + 1
				region_num += 1
				temp_array = DescribeLists(depthlist)
				averages = str(round(temp_array.average, 2))
				mediandepth = str(round(temp_array.median, 2))
				maxdepth = str(round(temp_array.max, 2))
				mindepth = str(round(temp_array.min, 2))
				beddepth.write("\t".join([chrom, rows[1], rows[2], averages, mediandepth, maxdepth, mindepth]) + '\n')
				for num in count_threshold:
					num = int(num)
					if temp_array.average >= num:
						regiondict[num] += 1
		uncover_range = list()
		for chrom in sort_bed_name_list(chroms.keys()):
			depthlist = chroms[chrom].values()
			temp_array = DescribeLists(depthlist)
			averages = str(round(temp_array.average, 2))
			mediandepth = str(round(temp_array.median, 2))
			maxdepth = str(round(temp_array.max, 2))
			mindepth = str(round(temp_array.min, 2))
			chromcover = [str(round(float(numpy.count_nonzero(temp_array.array >= i)) /
			                        len(temp_array.array) * 100, 2)) + "%" for i in count_threshold]
			chromdepth.write("\t".join([chrom, averages, mediandepth, maxdepth, mindepth] + chromcover) + '\n')
			if stepper == 'all':
				uncover_bases = [int(p) for p, d in chroms[chrom].iteritems() if d < uncover_threshold]
				uncover_range.extend(formact_number_list_to_range(uncover_bases, tag=chrom))
		if stepper == 'all':
			uncoverout = open(fileout + '.uncover.bed', 'w')
			uncoverout.writelines(uncover_range)
			uncoverout.close()
		coverout.write("Total Region Count : %s\n" % str(region_num))
		coverout.write("Len of region : %s\n" % str(total_base))
		for number in count_threshold:
			coverout.write("Coverage (>={0:s}X) : {1:s}%\n".format(
				str(number), str(round((float(rangedict[int(number)]) / region_num) * 100, 2))))
		for number in count_threshold:
			coverout.write("RegionCoverage (>={0:s}X) : {1:s}%\n".format(
				str(number), str(round((float(regiondict[int(number)]) / region_num) * 100, 2))))
		depths_out.close()
		coverout.close()
		beddepth.close()
		chromdepth.close()
		try:
			pysam.tabix_index(fileout, seq_col=0, start_col=1, end_col=1)
		except Exception as err:
			warnings.warn("fail to create index for %s : %s" % (fileout, err))


def sort_bed_name_list(bedlist, reverse=False):
	try:
		resultlist = sorted(bedlist, key=lambda s: int(s[3:]) if s[3:].isdigit() else sum([ord(n) for n in s[3:]]),
		                    reverse=reverse)
		return resultlist
	except Exception as err:
		raise Exception(err)


def formact_number_list_to_range(numberlits, tag=None):
	if not len(numberlits):
		return numberlits
	numberlits = sorted(numberlits)
	rangelist = list()
	for k, g in itertools.groupby(enumerate(numberlits), lambda (x, y): y - x):
		g = list(g)
		if tag is None:
			tag = str(k)
		final = "\t".join([tag, str(g[0][1]), str(g[-1][1])]) + '\n'
		rangelist.append(final)
	return rangelist


class DescribeLists(object):
	def __init__(self, numberlist):
		self.array = numpy.array(numberlist)
		self.average = numpy.mean(self.array)
		self.median = numpy.median(self.array)
		self.max = numpy.max(self.array)
		self.min = numpy.min(self.array)


def exon_plot(depthfile, trans, messages, outdir=os.getcwd(), prefix="Test", uncoverthreshold=4):
	depth = pysam.TabixFile(depthfile)
	genes = str
	plot_data = defaultdict(list)
	for exon_message in sorted(messages, key=lambda x: x[2]):
		if len(exon_message) != 7:
			continue
		genes, strand, exon_name, codingregion, refchrom, refstart, refend = exon_message
		if exon_name.startswith("IVS"):
			continue
		try:
			chrom = str(refchrom)
			start = int(refstart)
			end = int(refend)
		except Exception as err:
			warnings.warn("In line %s: %s" % ("\t".join(exon_message), err))
			continue
		depths = [int(depth_m.strip().split('\t')[2]) for depth_m in depth.fetch(chrom, start - 1, end)]
		if not len(depths):
			depths = [0] * (end - start + 1)
		depth_message = DescribeLists(depths)
		plot_data["averagedepth"].append(depth_message.average)
		plot_data["mediandepth"].append(depth_message.median)
		coverate = round((float(numpy.count_nonzero(depth_message.array >= uncoverthreshold)) / len(depths)) * 100, 2)
		plot_data["coverate"].append(coverate)
		plot_data["x_axis"].append(exon_name)
	depth.close()
	prefixname = ".".join([prefix, genes, trans])
	output = os.path.join(outdir, prefixname + '.png')
	size = max(len(plot_data["x_axis"]), 8)
	fig, axes = plt.subplots(2, 1)
	plt.suptitle(prefixname, fontsize=18)
	depth_ax, cover_ax = axes
	depth_ax.get_xaxis().set_major_locator(MaxNLocator(nbins=len(plot_data["x_axis"]) + 2, prune='upper'))
	depth_ax.bar(numpy.arange(len(plot_data["x_axis"])) - 0.2, numpy.array(plot_data["averagedepth"]), width=0.4,
	             label="Mean depth", align="center", color='c')
	depth_ax.bar(numpy.arange(len(plot_data["x_axis"])) + 0.2, numpy.array(plot_data["mediandepth"]), width=0.4,
	             color="pink", label="Median depth", align="center")
	depth_ax.set_title('Depths on all exon regions of %s' % trans, fontsize=15)
	depth_ax.set_ylabel("Depth (X)", fontsize=13)
	labels = numpy.array([''] + plot_data["x_axis"] + [''])
	depth_ax.set_xticklabels(labels, rotation=45, rotation_mode="anchor", size=min(180.0 / size, 13))
	depth_ax.legend(loc='best', prop={'size': 10})
	cover_ax.get_xaxis().set_major_locator(MaxNLocator(nbins=len(plot_data["x_axis"]) + 2, prune='upper'))
	cover_ax.plot(numpy.arange(len(plot_data["x_axis"])) + 0.2, numpy.array(plot_data["coverate"]),
	              'yo-', label="Coverage radio")
	cover_ax.set_title('Coverages on all exon regions of %s' % trans, fontsize=15)
	cover_ax.set_ylabel("Radio (%)", fontsize=13)
	cover_ax.set_ylim(0, 105)
	cover_ax.set_xticklabels(labels, rotation=45, rotation_mode="anchor", size=min(180.0 / size, 13))
	cover_ax.get_xaxis().set_visible(False)
	cover_ax.legend(loc='best', prop={'size': 10})
	plt.savefig(output, bbox_inches='tight', dpi=240, figsize=(size * 2, size))
	plt.close()


def graph_exon(depthfile, regiondict, outdir=os.getcwd(), prefix="Test"):
	import multiprocessing
	pool = multiprocessing.Pool(processes=min(multiprocessing.cpu_count(), len(regiondict.keys())))
	demo = "NM_000151.3"
	exon_plot(depthfile, demo, regiondict[demo], outdir, prefix)
	for transcript, trans_messages in regiondict.iteritems():
		pool.apply_async(exon_plot, (depthfile, demo, trans_messages, outdir, prefix))
	pool.close()
	pool.join()


def main():
	usage = "Usage: %prog [-h] [--version] --bam [bamfile] --config [config.json] -o [outdir]"
	author = 'JoeyHwong (joeyhwong@hotmail.com)'
	parser = OptionParser(usage=usage, version="%prog 1.0", epilog=author)
	parser.add_option('-o', '--outdir', help='The output dir [ default: %s ]' % os.getcwd(),
	                  dest='outdir', default=os.getcwd())
	parser.add_option('-b', '--bam', help='The bamfile (required) ', dest='bamfile')
	parser.add_option('-c', '--config', help='The config file (required) ', dest='config')
	parser.add_option('-p', '--plot', help='Plot exon depth graph', dest='plot', default=False, action="store_true")
	(options, args) = parser.parse_args()
	if not options.bamfile or not options.config:
		parser.print_help()
		return 1
	bamfile = os.path.realpath(options.bamfile)
	configfile = os.path.realpath(options.config)
	outdir = os.path.realpath(options.outdir)
	plot = options.plot
	baseconfig = ConfigPhrase(configfile)
	if baseconfig.config_stat == 0:
		raise ValueError('Sorry, expected arguments in config lost !!!')
	prepair_dir = CreateFolder(outdir)
	prepair_dir.create_subdirectory(["bed_anno", "depth_anno", "depth_anno/Target", "depth_anno/Flank"])
	tools = baseconfig.tools
	bedanno = BedAnno(tools)
	bam = BamDepth(bamfile)
	sample_name = bam.sample_name
	with open(os.path.join(outdir, 'bed_anno', sample_name + '.Target.bed.anno'), 'w') as f_out:
		f_out.writelines(bedanno.target_anno)
	with open(os.path.join(outdir, 'bed_anno', sample_name + '.Flank.bed.anno'), 'w') as f_out:
		f_out.writelines(bedanno.flank_anno)

	for bedtype in ('Target', 'Flank'):
		bam.countdepthcover(tools["%s bed" % bedtype],
		                    os.path.join(outdir, 'depth_anno', bedtype, sample_name + '.RmdupDepth'))
		bam.countdepthcover(tools["%s bed" % bedtype],
		                    os.path.join(outdir, 'depth_anno', bedtype, sample_name + '.RawDepth'),
		                    'nofilter')
	uncoverbed = glob.glob(os.path.join(outdir, 'depth_anno', "*", sample_name + ".RmdupDepth.uncover.bed"))
	for bed in uncoverbed:
		uncover_out = re.sub("\.bed$", ".anno", bed)
		bedanno.uncoveranno(bed, uncover_out)
	if plot:
		prepair_dir.create_subdirectory(["graph_exon"])
		graph_exon(os.path.join(outdir, 'depth_anno', 'Flank', sample_name + '.RmdupDepth.gz'),
		           bedanno.trans_region, os.path.join(outdir, 'graph_exon'), sample_name)


if __name__ == '__main__':
	sys.exit(main())
