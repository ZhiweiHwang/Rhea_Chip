#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
Phasing the vcf file, add NB, PB
"""

import pysam
import os
import re
import sys
import gzip
import itertools
import tempfile
import shutil


def split_format(string1, string2):
	fmat, val = string1.split(":"), string2.split(":")
	val = val[:len(fmat)] if len(fmat) < len(val) else val + ['.'] * (len(fmat) - len(val))
	return fmat, val


def add_marker(string1, string2, key1, key2):
	fmat, val = split_format(string1, string2)
	string2 = ":".join(val)
	if key1 in fmat:
		return string1, string2
	else:
		return string1 + ':' + key1, string2 + ':' + key2


def support_reads(bamfile, chrom, pos, ref, alt, gtype):
	pybamfile = pysam.AlignmentFile(bamfile, 'rb')
	pos = int(pos)
	return_reads = {'first': list(), 'last': list()}
	if len(ref) > len(alt):
		ref = '-' * (len(ref) - len(alt))
	g1, tag, g2 = re.compile('(\d)([/\|])(\d)').match(gtype).groups()
	if g1 == 1:
		refer = 'last'
		mut = 'first'
	else:
		refer = 'first'
		mut = 'last'
	for pc in pybamfile.pileup(chrom, pos - 1, pos):
		if pc.pos != pos - 1:
			continue
		for pr in pc.pileups:
			align = pr.alignment
			if pr.is_del or pr.indel < 0:
				bases = '-' * abs(pr.indel)
			else:
				pos_start = pr.query_position
				pos_end = pr.query_position + pr.indel + 1
				bases = align.seq[pos_start:pos_end]
			if bases == ref:
				return_reads[refer].append(align.query_name)
			elif bases == alt:
				return_reads[mut].append(align.query_name)
	return return_reads


def reads_phasing(reads1, reads2):
	compare1 = len([reads for reads in reads1['first'] if reads in reads2['first']])
	compare2 = len([reads for reads in reads1['first'] if reads in reads2['last']])
	compare3 = len([reads for reads in reads1['last'] if reads in reads2['first']])
	compare4 = len([reads for reads in reads1['last'] if reads in reads2['last']])
	if compare1 > compare2 and compare4 > compare3:
		return 1
	if compare2 > compare1 and compare3 > compare4:
		return 2
	return 0


def readfile(infile, bamfile):
	returnlist = list()
	addformat = "##FORMAT=<ID=NB,Number=1,Type=String,Description=\"Neighbor Block\">\n" \
	            "##FORMAT=<ID=PB,Number=1,Type=String,Description=\"Phasing Block\">\n"
	f_in = gzip.open(infile, 'rb').readlines() if infile.endswith('.gz') else open(infile, 'rb').readlines()
	vcfdict = {k: list(g) for k, g in itertools.groupby(f_in, lambda x: x.startswith('#'))}
	header, variantinfo = vcfdict[True], [line.strip().split('\t') for line in vcfdict[False]]
	formatline = max([num for num, line in enumerate(header) if line.startswith('##FORMAT')])
	header.insert(formatline + 1, addformat)
	nbnumber = {str(i): 0 for i in range(1, 26)}
	pbnumber = {str(i): 0 for i in range(1, 26)}
	pre_type, rev, pre_haplo = None, 0, dict()
	for num, rows in enumerate(variantinfo[1:], 1):
		if len(rows) != 10:
			continue
		pos = int(rows[1]) + len(rows[3]) - 1
		chrom = rows[0].lower().replace('chr', '').replace('x', '23').replace('y', '24')
		chrom = re.sub(".*m.*", "25", chrom)
		mut = min(rows[4].split(","), key=len)
		ref = rows[3]
		if rows[0].lower() != variantinfo[num - 1][0].lower():
			pre_type, rev, pre_haplo = None, 0, dict()
		muttypes = 'SNP' if len(ref) == 1 and len(mut) == 1 else 'INDEL'
		formats, values = split_format(rows[-2], rows[-1])
		preformat, prevalue = split_format(variantinfo[num - 1][-2], variantinfo[num - 1][-1])
		formatdict = dict(zip(formats, values))
		preformatdict = dict(zip(preformat, prevalue))
		gtype = re.compile("(\d)([/\|])(\d)").match(formatdict['GT'])
		pregtype = re.compile("(\d)([/\|])(\d)").match(preformatdict['GT'])
		if gtype and (rows[6] == 'PASS' or rows[6] == '.'):
			g1, tag, g2 = gtype.groups()
			if g1 == g2:
				formatdict['GT'] = str(g1) + '|' + str(g2)
			elif 'HP' in formatdict:
				if 'HP' in preformatdict:
					prehaplo = re.split(',|-', preformatdict['HP'])
					preg1, pretag, preg2 = pregtype.groups()
					if prehaplo[0] not in pre_haplo:
						pre_haplo[prehaplo[0]] = dict()
					if prehaplo[1] not in pre_haplo[prehaplo[0]]:
						pre_haplo[prehaplo[0]][prehaplo[1]] = preg1
					if prehaplo[2] not in pre_haplo:
						pre_haplo[prehaplo[2]] = dict()
					if prehaplo[3] not in pre_haplo[prehaplo[2]]:
						pre_haplo[prehaplo[2]][prehaplo[3]] = preg2
				haplo = re.split(',|-', formatdict['HP'])
				if haplo[0] in pre_haplo and haplo[1] in pre_haplo[haplo[0]] and haplo[3] in pre_haplo[haplo[0]]:
					formatdict['GT'] = str(pre_haplo[haplo[0]][haplo[1]]) + '|' + str(pre_haplo[haplo[0]][haplo[3]])
					if ':PB' not in variantinfo[num - 1][-2]:
						pbnumber[chrom] += 1
					variantinfo[num][-2], variantinfo[num][-1] = add_marker(variantinfo[num][-2],
					                                                        variantinfo[num][-1],
					                                                        'PB',
					                                                        "P{:0>2d}{:0>5d}".format(
						                                                        int(chrom), int(pbnumber[chrom])))
					variantinfo[num - 1][-2], variantinfo[num - 1][-1] = add_marker(variantinfo[num - 1][-2],
					                                                                variantinfo[num - 1][-1],
					                                                                'PB',
					                                                                "P{:0>2d}{:0>5d}".format(
						                                                                int(chrom),
						                                                                int(pbnumber[chrom])))
				else:
					pre_haplo[haplo[0]] = dict()
					pre_haplo[haplo[0]][haplo[1]] = str(g1)
					pre_haplo[haplo[0]][haplo[3]] = str(g2)
			elif muttypes == 'SNP':
				if pre_type == 'SNP':
					if g1 != g2:
						if tag == '|':
							if ':PB' not in variantinfo[num - 1][-2]:
								pbnumber[chrom] += 1
							variantinfo[num][-2], variantinfo[num][-1] = add_marker(variantinfo[num][-2],
							                                                        variantinfo[num][-1],
							                                                        'PB',
							                                                        "P{:0>2d}{:0>5d}".format(
								                                                        int(chrom),
								                                                        int(pbnumber[chrom])))
							variantinfo[num - 1][-2], variantinfo[num - 1][-1] = add_marker(variantinfo[num - 1][-2],
							                                                                variantinfo[num - 1][-1],
							                                                                'PB',
							                                                                "P{:0>2d}{:0>5d}".format(
								                                                                int(chrom),
								                                                                int(pbnumber[chrom])))
							if rev:
								formatdict['GT'] = str(g2) + '|' + str(g1)
						else:
							rev = 0
				elif pre_type == 'INDEL':
					if g1 != g2:
						pre_reads = support_reads(bamfile, variantinfo[num - 1][0], variantinfo[num - 1][1],
						                          variantinfo[num - 1][3], variantinfo[num - 1][4], preformatdict['GT'])
						now_reads = support_reads(bamfile, variantinfo[num][0], variantinfo[num][1],
						                          variantinfo[num][3], variantinfo[num][4], formatdict['GT'])
						phasing = reads_phasing(pre_reads, now_reads)
						if phasing:
							if ':PB' not in variantinfo[num - 1][-2]:
								pbnumber[chrom] += 1
							variantinfo[num][-2], variantinfo[num][-1] = add_marker(variantinfo[num][-2],
							                                                        variantinfo[num][-1], 'PB',
							                                                        "P{:0>2d}{:0>5d}".format(
								                                                        int(chrom),
							                                                            int(pbnumber[chrom])))
							variantinfo[num - 1][-2], variantinfo[num - 1][-1] = add_marker(variantinfo[num - 1][-2],
							                                                                variantinfo[num - 1][-1],
							                                                                'PB',
							                                                                "P{:0>2d}{:0>5d}".format(
								                                                                int(chrom),
								                                                                int(pbnumber[chrom])))
							if phasing == 1:
								rev = 0 if rev == 0 else 1
							else:
								rev = 1 if rev == 0 else 0
							formatdict['GT'] = str(g1) + '|' + str(g2) if rev == 0 else str(g2) + '|' + str(g1)
						else:
							rev = 0
			else:
				if g1 != g2:
					pre_reads = support_reads(bamfile, variantinfo[num - 1][0], variantinfo[num - 1][1],
					                          variantinfo[num - 1][3], variantinfo[num - 1][4], preformatdict['GT'])
					now_reads = support_reads(bamfile, variantinfo[num][0], variantinfo[num][1],
					                          variantinfo[num][3], variantinfo[num][4], formatdict['GT'])
					phasing = reads_phasing(pre_reads, now_reads)
					if phasing:
						if ':PB' not in variantinfo[num - 1][-2]:
							pbnumber[chrom] += 1
						variantinfo[num][-2], variantinfo[num][-1] = add_marker(variantinfo[num][-2],
						                                                        variantinfo[num][-1], 'PB',
						                                                        "P{:0>2d}{:0>5d}".format(int(chrom),
						                                                                                 int(pbnumber[
							                                                                                 chrom])))
						variantinfo[num - 1][-2], variantinfo[num - 1][-1] = add_marker(variantinfo[num - 1][-2],
						                                                                variantinfo[num - 1][-1], 'PB',
						                                                                "P{:0>2d}{:0>5d}".format(
							                                                                int(chrom),
							                                                                int(pbnumber[chrom])))
						if phasing == 1:
							rev = 0 if rev == 0 else 1
						else:
							rev = 1 if rev == 0 else 0
						formatdict['GT'] = str(g1) + '|' + str(g2) if rev == 0 else str(g2) + tag + str(g1)
					else:
						rev = 0
			pre_type = muttypes

		if rows[0] == variantinfo[num - 1][0] and \
			pos - (int(variantinfo[num - 1][1]) + len(variantinfo[num - 1][3]) - 1) <= 5:
			if ':NB' not in variantinfo[num - 1][-2]:
				nbnumber[chrom] += 1
			variantinfo[num][-2], variantinfo[num][-1] = add_marker(variantinfo[num][-2], variantinfo[num][-1], 'NB',
			                                                        "N{:0>2d}{:0>5d}".format(int(chrom),
			                                                                                 int(nbnumber[chrom])))
			variantinfo[num - 1][-2], variantinfo[num - 1][-1] = add_marker(variantinfo[num - 1][-2],
			                                                                variantinfo[num - 1][-1], 'NB',
			                                                                "N{:0>2d}{:0>5d}".format(
				                                                                int(chrom), int(nbnumber[chrom])))
		formatlast, valuelast = split_format(variantinfo[num][-2], variantinfo[num][-1])
		for i in range(len(formatlast)):
			if formatlast[i] in formatdict:
				valuelast[i] = formatdict[formatlast[i]]
		variantinfo[num][-2], variantinfo[num][-1] = ":".join(formatlast), ":".join(valuelast)
	returnlist.extend(header)
	returnlist.extend(["\t".join(i) + '\n' for i in variantinfo])
	return returnlist


def main():
	usage = 'python %s <bamfile> <VcfIn> <VcfOut>' % sys.argv[0]
	if len(sys.argv) != 4:
		return usage
	bamfile = os.path.realpath(sys.argv[1])
	infile = os.path.realpath(sys.argv[2])
	outfile = os.path.realpath(sys.argv[3])
	temppath = tempfile.mkdtemp()
	resultfileout = os.path.join(temppath, os.path.basename(outfile))
	f_out = open(re.sub("\.gz$", "", resultfileout), 'w')
	f_out.writelines(readfile(infile, bamfile))
	f_out.close()
	if resultfileout.endswith('.gz'):
		try:
			resultfileout = pysam.tabix_index(re.sub("\.gz$", "", resultfileout), preset='vcf', force=True)
			shutil.copy(resultfileout + '.tbi', outfile + '.tbi')
		except Exception:
			pass
	shutil.copy(resultfileout, outfile)

if __name__ == '__main__':
	sys.exit(main())
