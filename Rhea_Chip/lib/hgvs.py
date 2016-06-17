#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
import pysam
import re
import sqlite3
import sys
from smart_open import smart_open
from collections import defaultdict, namedtuple


def get_pair(seq):
	rule = {'A': 'T', 'C': 'G', 'B': 'V', 'D': 'H', 'G': 'C', 'H': 'D', 'K': 'M', 'M': 'K', 'S': 'S', 'R': 'Y',
	        'T': 'A', 'W': 'W', 'V': 'B', 'Y': 'R', 'a': 't', 'c': 'g', 'b': 'v', 'd': 'h', 'g': 'c', 'h': 'd',
	        'k': 'm', 'm': 'k', 's': 's', 'r': 'y', 't': 'a', 'w': 'w', 'v': 'b', 'y': 'r'}
	patten = re.compile("([^ATCGatcgRYKMSWBDHVrykmswbdhv,]+)?([ATCGatcgRYKMSWBDHVrykmswbdhv]+)"
	                    "([^,ATCGatcgRYKMSWBDHVrykmswbdhv]+)?").findall(seq)
	if not patten:
		return seq
	new = list()
	for s in patten:
		new.append("{0}{1}{2}".format(s[0], "".join(map(lambda x: rule[x], s[1]))[::-1], s[2]))
	return ",".join(new)


class HGVSName(object):
	pass


class HGVSRegex(object):
	"""
	All regular expression for HGVS names.
	"""

	# DNA syntax
	# http://www.hgvs.org/mutnomen/standards.html#nucleotide
	BASE = "[ATCGatcgnN]|\d+"
	BASES = "[ATCGatcgnN]+|\d+"
	DNA_REF = "(?P<ref>" + BASES + ")"
	DNA_ALT = "(?P<alt>" + BASES + ")"

	# Mutation types
	EQUAL = "(?P<mutation_type>=)"
	SUB = "(?P<mutation_type>>)"
	INS = "(?P<mutation_type>ins)"
	DEL = "(?P<mutation_type>del)"
	DUP = "(?P<mutation_type>dup)"

	# Simple coordinate syntax
	COORD_START = "(?P<start_hgvs>\d+)"
	COORD_END = "(?P<end_hgvs>\d+)"
	COORD_RANGE = COORD_START + "_" + COORD_END

	# cDNA coordinate syntax
	CDNA_START = "(?P<start_hgvs>[-\*]?\d+)(?P<start_offset_prefix>-|\+)?(?P<start_offset>\d+)?"
	CDNA_END = r"(?P<end_hgvs>[-\*]?\d+)(?P<end_offset_prefix>-|\+)?(?P<end_offset>\d+)?"
	CDNA_RANGE = CDNA_START + "_" + CDNA_END

	# cDNA allele syntax
	CDNA_ALLELE = [  # No change
		CDNA_START + DNA_REF + EQUAL,  # Substitution
		CDNA_START + DNA_REF + SUB + DNA_ALT,  # 1bp insertion, deletion, duplication
		CDNA_START + INS + DNA_ALT,
		CDNA_START + DEL + DNA_REF,
		CDNA_START + DUP + DNA_REF,
		CDNA_START + DEL,
		CDNA_START + DUP,  # Insertion, deletion, duplication
		CDNA_RANGE + INS + DNA_ALT,
		CDNA_RANGE + DEL + DNA_REF,
		CDNA_RANGE + DUP + DNA_REF,
		CDNA_RANGE + DEL,
		CDNA_RANGE + DUP,  # Indels
		"(?P<delins>" + CDNA_START + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
		"(?P<delins>" + CDNA_RANGE + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
		"(?P<delins>" + CDNA_START + 'delins' + DNA_ALT + ")",
		"(?P<delins>" + CDNA_RANGE + 'delins' + DNA_ALT + ")",
	]

	CDNA_ALLELE_REGEXES = [re.compile("^" + regex + "$")
	                       for regex in CDNA_ALLELE]

	# IVS allele syntax
	IVS_COORD = "([ad]s)?(?P<offset_prefix>-|\+)(?P<offset>\d+)"
	IVS_ALLELE = [
		COORD_START + IVS_COORD + DNA_REF + EQUAL,
		COORD_START + IVS_COORD + DNA_REF + SUB + DNA_ALT,
		COORD_START + IVS_COORD + INS + DNA_ALT,
		COORD_START + IVS_COORD + DEL + DNA_REF,
		COORD_START + IVS_COORD + DUP + DNA_REF,
		COORD_START + IVS_COORD + DEL,
		COORD_START + IVS_COORD + DUP,
		"(?P<delins>" + COORD_START + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
		"(?P<delins>" + COORD_START + 'delins' + DNA_ALT + ")",
	]
	IVS_ALLELE_REGEXES = [re.compile("^" + regex + "$")
	                      for regex in IVS_ALLELE]

	# Genomic allele syntax
	GENOMIC_ALLELE = [  # No change
		COORD_START + DNA_REF + EQUAL,  # Substitution
		COORD_START + DNA_REF + SUB + DNA_ALT,  # 1bp insertion, deletion, duplication
		COORD_START + INS + DNA_ALT,
		COORD_START + DEL + DNA_REF,
		COORD_START + DUP + DNA_REF,
		COORD_START + DEL,
		COORD_START + DUP,  # Insertion, deletion, duplication
		COORD_RANGE + INS + DNA_ALT,
		COORD_RANGE + DEL + DNA_REF,
		COORD_RANGE + DUP + DNA_REF,
		COORD_RANGE + DEL,
		COORD_RANGE + DUP,  # Indels
		"(?P<delins>" + COORD_START + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
		"(?P<delins>" + COORD_RANGE + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
		"(?P<delins>" + COORD_START + 'delins' + DNA_ALT + ")",
		"(?P<delins>" + COORD_RANGE + 'delins' + DNA_ALT + ")",
	]

	GENOMIC_ALLELE_REGEXES = [re.compile("^" + regex + "$")
	                          for regex in GENOMIC_ALLELE]

	def parse_hgvs(self, chgvs):
		hgvs = HGVSName()
		chgvs = re.sub("\s+", "", chgvs.lower())
		if chgvs.startswith("ivs"):
			hgvs.kind = "IVS"
			details = chgvs[3:]
			for regex in self.IVS_ALLELE_REGEXES:
				match = re.match(regex, details)
				if match:
					groups = match.groupdict()
					if groups.get('delins'):
						hgvs.mutation_type = 'delins'
					else:
						hgvs.mutation_type = groups['mutation_type']
					hgvs.cdna_start = groups.get('start_hgvs')
					hgvs.offset_prefix = groups.get('offset_prefix')
					hgvs.offset = groups.get('offset')
					hgvs.ref_allele = groups.get('ref', '').upper()
					hgvs.alt_allele = groups.get('alt', '').upper()
					if hgvs.ref_allele.isdigit():
						hgvs.ref_allele = "N" * int(hgvs.ref_allele)
					if hgvs.alt_allele.isdigit():
						hgvs.alt_allele = "N" * int(hgvs.alt_allele)
					if hgvs.mutation_type == "=":
						hgvs.alt_allele = hgvs.ref_allele
					return hgvs
		elif chgvs.startswith("g."):
			hgvs.kind = "genomic"
			details = chgvs[2:]
			for regex in self.GENOMIC_ALLELE_REGEXES:
				match = re.match(regex, details)
				if match:
					groups = match.groupdict()
					if groups.get('delins'):
						hgvs.mutation_type = 'delins'
					else:
						hgvs.mutation_type = groups['mutation_type']
					hgvs.start = int(groups.get('start_hgvs'))
					hgvs.end = int(groups.get('end_hgvs')) if groups.get('end_hgvs') else hgvs.start
					hgvs.ref_allele = groups.get('ref', '').upper()
					hgvs.alt_allele = groups.get('alt', '').upper()
					if hgvs.ref_allele.isdigit():
						hgvs.ref_allele = "N" * int(hgvs.ref_allele)
					if hgvs.alt_allele.isdigit():
						hgvs.alt_allele = "N" * int(hgvs.alt_allele)
					if hgvs.mutation_type == "=":
						hgvs.alt_allele = hgvs.ref_allele
					return hgvs
		else:
			if re.compile("(-)?\d+").match(chgvs):
				chgvs = "c." + str(chgvs)
			hgvs.kind = "CDNA"
			if chgvs.startswith("m."):
				hgvs.kind = "MT"
			details = chgvs[2:]
			for regex in self.CDNA_ALLELE_REGEXES:
				match = re.match(regex, details)
				if match:
					groups = match.groupdict()
					if groups.get('delins'):
						hgvs.mutation_type = 'delins'
					else:
						hgvs.mutation_type = groups['mutation_type']
					hgvs.cdna_start = groups.get('start_hgvs')
					hgvs.cdna_start_offset_prefix = groups.get('start_offset_prefix')
					hgvs.cdna_start_offset = groups.get('start_offset')
					if 'end_hgvs' not in groups:
						hgvs.cdna_end = hgvs.cdna_start
						hgvs.cdna_end_offset_prefix = hgvs.cdna_start_offset_prefix
						hgvs.cdna_end_offset = hgvs.cdna_start_offset
					else:
						hgvs.cdna_end = groups.get('end_hgvs')
						hgvs.cdna_end_offset_prefix = groups.get('end_offset_prefix')
						hgvs.cdna_end_offset = groups.get('end_offset')
					hgvs.ref_allele = groups.get('ref', '').upper()
					hgvs.alt_allele = groups.get('alt', '').upper()
					if hgvs.ref_allele.isdigit():
						hgvs.ref_allele = "N" * int(hgvs.ref_allele)
					if hgvs.alt_allele.isdigit():
						hgvs.alt_allele = "N" * int(hgvs.alt_allele)
					if hgvs.mutation_type == "=":
						hgvs.alt_allele = hgvs.ref_allele
					return hgvs


class HGVS(object):
	def __init__(self, refgene, reference=None, trans=None, genes=None):
		self.db = os.path.abspath(refgene)
		self.reference = os.path.abspath(reference) if reference else None
		if not os.path.isfile(self.db):
			raise IOError("Expected db lost !!!")
		self.dbref = sqlite3.connect(self.db)
		self.refer = pysam.FastaFile(self.reference) if self.reference else None
		self.trans = self.get_primary(trans, ext=True)
		self.genes = self.get_primary(genes)
		self.hgvsname = namedtuple("HGVSMessage", ("Transcript", "cHgvs", "PrimaryTag", "geneSym", "Protein",
		                                           "Stand", "blockAttr", "ExonRegions", "Chrom", "Start",
		                                           "Stop", "Refer", "Call"))

	def __del__(self):
		self.dbref.close()
		if self.refer is not None:
			self.refer.close()

	@staticmethod
	def get_primary(f_in, ext=False):
		if f_in is not None:
			f_in = os.path.abspath(f_in)
		else:
			return None
		subset = set()
		with smart_open(f_in) as f:
			for line in f:
				line = line.strip().split()[-1]
				if ext:
					line = re.sub("\.[^\.]+$", "", line)
				subset.add(line)
		return subset

	def split_regions(self, chrom, start, end):
		if chrom.startswith("chrM"):
			chrom = 'chrM_NC_012920.1'
		feback = defaultdict(set)
		start = int(start)
		end = int(end)
		sql = 'SELECT primaryTag, Trans, Version, geneSym, txStart, txEnd FROM RefHeader WHERE Chrom="{0}" and ' \
		      '((txStart<={1} and txEnd>={1}) or (txStart>={1} and txEnd<={2}) or (txStart<={2} and txEnd>={2}))'
		sencur = self.dbref.execute(sql.format(chrom, start, end)).fetchall()
		if sencur is None:
			feback[(chrom, start, end)] = 'Intergenic'
			return feback
		choose_tran = dict()
		nor_len = defaultdict(int)
		tran_dict = dict()
		trans_list = set()
		for primay, tran, version, gene, tx_start, tx_end in sencur:
			if self.trans is not None and tran not in self.trans:
				continue
			if self.genes is not None and gene not in self.genes:
				continue
			tx_start = int(tx_start)
			tx_end = int(tx_end)
			length = tx_end - tx_start + 1
			tran = ".".join([tran, version])
			if self.trans is None:
				trans_list.add(tran)
			else:
				trans_list = {tran}
			if chrom.startswith('chrM'):
				gene = "MT-" + re.sub("^MT-", "", gene)
			keys = gene if chrom.startswith("chrM") else tran
			tran_dict[tran] = (keys, tx_start, tx_end)
			if length > nor_len[gene]:
				choose_tran[gene] = tran
				nor_len[gene] = length
			if primay == 'Y':
				choose_tran[gene] = tran
				nor_len[gene] = sys.maxint
		final = sorted([tran_dict[tran] for gene, tran in choose_tran.iteritems()], key=lambda x: x[1])
		for tran, begin, stop in final:
			if begin > start:
				feback[(chrom, start, begin - 1)] = 'Intergenic'
			else:
				begin = start
			stop = min(stop, end)
			if begin > stop:
				continue
			feback[(chrom, begin, stop)] = tran
			start = stop + 1
		if start < end:
			feback[(chrom, start, end)] = 'Intergenic'
		return feback, trans_list

	def format_hgvs_name(self, chrom, pos, transcript, baselevel=0):
		refer_chr = re.sub("^chr", "", chrom)
		try:
			refer_chr = "NC_{:0>6}".format(int(refer_chr))
		except ValueError:
			if re.compile("m", re.I).match(refer_chr):
				refer_chr = 'NC_012920.1'
			elif re.compile("x", re.I).match(refer_chr):
				refer_chr = "NC_000023"
			elif re.compile("y", re.I).match(refer_chr):
				refer_chr = "NC_000024"
			else:
				refer_chr = chrom
		tmp_hgvs = (transcript, refer_chr, '.', 'Y', '+', 'g.%d' % pos, '.', '.')
		pos = int(pos) - baselevel + 1
		try:
			phrase_trans, phrase_version = transcript.split(".")
		except ValueError:
			phrase_trans, phrase_version = transcript, None
		if transcript.startswith("NM_"):
			cdsopt = 'c.'
		elif transcript.startswith("NR_"):
			cdsopt = 'n.'
		elif transcript.startswith("MT_"):
			cdsopt = 'm.'
		else:
			return tmp_hgvs
		selectsql = "SELECT rv.geneSym, rv.primaryTag, rv.Strand, rv.Protein, ms.blockAttr, ms.Regions, " \
		            "ms.chrom, ms.Start, ms.End, ms.hgvsStart, ms.hgvsEnd FROM Refcontent AS ms, " \
		            "RefHeader AS rv WHERE ms.chrom=\"{0}\" AND ms.Start<={1} AND ms.End >={1} AND ms.RefIds " \
		            "= rv.RefIds and rv.Trans = \"{2}\"".format(chrom, pos, phrase_trans)
		if phrase_version is not None:
			selectsql += "and rv.Version = \"%s\"" % phrase_version
		sencur = self.dbref.execute(selectsql)
		senres = sencur.fetchall()
		if senres is None or not len(senres):
			return None
		for rows in sorted(senres, key=lambda s: s[1]):
			gene, primary, strand, protein, block, funregion, lcs_chr, lcs_start, lcs_end, t_start, t_end = rows
			offset = pos - int(lcs_start)
			pre = "*" if re.compile("\*").match(t_start) else ''
			if pre:
				t_start = t_start[1:]
			if funregion.startswith("EX"):
				chgvs = int(t_start) + offset if strand == '+' else int(t_start) - offset
				return gene, transcript, protein, primary, strand, "".join([cdsopt, pre, str(chgvs)]), block, funregion
			elif funregion.startswith("IVS") or funregion == "promoter":
				if int(lcs_end) - pos > offset:
					temp = re.compile("(-?[0-9]+)([+-])?([0-9]+)?").match(t_start)
					if temp:
						t_start, l_strand, l_off = temp.groups()
						if l_off is not None:
							offset = int(l_off) + offset if l_strand == strand else int(l_off) - offset
							strand = l_strand
						if funregion == "promoter":
							promoter_hgvs = int(t_start) + offset if strand == '+' else int(t_start) - offset
							chgvs = "".join([cdsopt, pre, str(promoter_hgvs)])
						else:
							chgvs = "".join([cdsopt, pre, str(t_start), strand, str(offset)])
					else:
						continue
				else:
					offset = int(lcs_end) - pos
					temp = re.compile("(-?[0-9]+)([+-])?([0-9]+)?").match(t_end)
					if temp:
						t_end, l_strand, l_off = temp.groups()
						if l_off is not None:
							offset = int(l_off) - offset if l_strand == strand else int(l_off) + offset
						else:
							offset += 1
							l_strand = '-' if strand == "+" else "+"
						if funregion == "promoter":
							promoter_hgvs = int(t_end) + offset if l_strand == '+' else int(t_end) - offset
							chgvs = "".join([cdsopt, pre, str(promoter_hgvs)])
						else:
							chgvs = "".join([cdsopt, pre, str(t_end), l_strand, str(offset)])
					else:
						continue
				return gene, transcript, protein, primary, strand, chgvs, block, funregion
			else:
				continue
		return None

	def mutalyzer(self, chrom, start, stop, refer=None, alter=None, tran_list=None):
		feback = set()
		bed_anno = True if refer is None and alter is None else False
		if refer is None:
			refer = str()
		if alter is None:
			alter = str()
		for tran in tran_list:
			s_hgvs = self.format_hgvs_name(chrom, start, tran)
			if s_hgvs is None:
				continue
			s_g, s_t, s_p, s_pr, s_s, s_c, s_b, s_f = s_hgvs
			_, _, _, _, _, e_c, e_b, e_f = s_hgvs
			prefix_hgvs = s_c
			if s_s == "-":
				refer = get_pair(refer)
				alter = get_pair(alter)
			if refer == alter and not bed_anno:
				chgvs = "%s%s=" % (prefix_hgvs, refer)
			elif len(refer) == len(alter) == 1:
				chgvs = "%s%s>%s" % (prefix_hgvs, refer, alter)
			else:
				if stop - 1 != start:
					e_hgvs = self.format_hgvs_name(chrom, stop, tran, 1)
					if e_hgvs is None:
						continue
					_, _, _, _, _, e_c, e_b, e_f = e_hgvs
					if start >= stop:
						s_c, e_c = e_c, s_c
					if s_s == '+':
						s_b = "{0}-{1}".format(s_b, e_b) if s_b != e_b else s_b
						s_f = "{0}-{1}".format(s_f, e_f) if s_f != e_f else s_f
					else:
						s_b = "{0}-{1}".format(e_b, s_b) if s_b != e_b else e_b
						s_f = "{0}-{1}".format(e_f, s_f) if s_f != e_f else e_f
					prefix_hgvs = "_".join([s_c, e_c[2:]]) if s_s == "+" else "_".join([e_c, s_c[2:]])
				if len(refer) == 0:
					chgvs = "%sins%s" % (prefix_hgvs, alter)
				elif len(alter) == 0:
					chgvs = "%sdel%s" % (prefix_hgvs, refer)
				else:
					chgvs = "%sdel%sins%s" % (prefix_hgvs, refer, alter)
			f = prefix_hgvs if bed_anno else chgvs
			final = (s_t, f, s_pr, s_g, s_p, s_s, s_b, s_f, chrom, start, stop, refer, alter)
			feback.add(self.hgvsname._make(final))
		return feback

	def annovar(self, chrom, pos, refer, alter):
		chrom = "chr" + re.sub("^chr", "", chrom)
		start = int(pos) - 1
		refer = refer.upper()
		alter = alter.upper()
		tmp_ref = refer
		while len(refer) and len(alter) and refer[0] == alter[0]:
			refer = refer[1:]
			alter = alter[1:]
			start += 1
		stop = start + len(refer)
		if not len(refer) and not len(alter):
			refer = alter = tmp_ref[0]
		anno_set = set()
		regions, tran_list = self.split_regions(chrom, start, stop)
		regions = sorted(regions.iteritems(), key=lambda s: int(s[0][1]))
		first_region, first_name = regions.pop()
		c, b, e = first_region
		r = refer[b - start: e - start]
		a = alter
		if first_name == 'Intergenic':
			anno_set |= self.mutalyzer(c, int(b), int(e), r, a, ['Intergenic'])
		else:
			anno_set |= self.mutalyzer(c, int(b), int(e), r, a, tran_list)
		for after_region, after_name in regions:
			c, b, e = after_region
			r = refer[b - start: e - start]
			a = str
			if after_name == 'Intergenic':
				anno_set |= self.mutalyzer(c, int(b), int(e), r, a, ['Intergenic'])
			else:
				anno_set |= self.mutalyzer(c, int(b), int(e), r, a, tran_list)
		return anno_set

	def annobed(self, chrom, start, stop):
		chrom = "chr" + re.sub("^chr", "", chrom)
		start = int(start) - 1
		stop = int(stop)
		anno_set = set()
		regions, tran_list = self.split_regions(chrom, start, stop)
		regions = sorted(regions.iteritems(), key=lambda s: int(s[0][1]))
		first_region, first_name = regions.pop()
		c, b, e = first_region
		if first_name == 'Intergenic':
			anno_set |= self.mutalyzer(c, int(b), int(e), None, None, ['Intergenic'])
		else:
			anno_set |= self.mutalyzer(c, int(b), int(e), None, None, tran_list)
		for after_region, after_name in regions:
			c, b, e = after_region
			if after_name == 'Intergenic':
				anno_set |= self.mutalyzer(c, int(b), int(e), None, None, ['Intergenic'])
			else:
				anno_set |= self.mutalyzer(c, int(b), int(e), None, None, tran_list)
		return anno_set

	def parse_hgvs_name(self, chgvs, transcript=None, genesymbol=None):
		assert transcript or genesymbol, "Trans or Gene must be set one"
		trans_set = set()
		if transcript:
			trans = transcript.split(".")
			trans_select = self.dbref.execute("SELECT * FROM RefHeader where Trans = \"%s\"" % trans[0])
			trans_set = set(trans_select.fetchall())
			if len(trans) > 1:
				new_set = filter(lambda s: s[2] != trans[1], trans_set)
				if len(new_set):
					trans_set = new_set
		if genesymbol and not len(trans_set):
			trans_select = self.dbref.execute("SELECT * FROM RefHeader where geneSym = \"%s\"" % genesymbol)
			trans_set = set(trans_select.fetchall())
		regex = HGVSRegex()
		hgvs = regex.parse_hgvs(chgvs)
		if hgvs is None:
			raise AttributeError("Input Formact Error: %s" % chgvs)
		feback = set()
		for refseq_rows in trans_set:
			ids, genesym, geneid, trans, t_version, primarytag, protein, strand, chrom = refseq_rows[:9]
			refseq_select = self.dbref.execute("SELECT Start, End, Regions, blockAttr, hgvsStart, hgvsEnd "
			                                   "FROM Refcontent WHERE RefIds=%d ORDER BY Start" % ids).fetchall()
			trans_name = ".".join([trans, t_version]) if hgvs.kind != "MT" else trans
			hgvs.chrom = chrom
			hgvs.geneSym = genesym
			hgvs.Entrez_id = geneid
			hgvs.PrimaryTag = primarytag
			hgvs.protein = protein
			hgvs.Strand = strand
			if hgvs.Strand == '-':
				hgvs.ref_allele = get_pair(hgvs.ref_allele)
				hgvs.alt_allele = get_pair(hgvs.alt_allele)
			if hgvs.kind == "IVS":
				ivs_number = "IVS" + hgvs.cdna_start
				max_p = -1
				min_p = sys.maxint
				for rows in refseq_select:
					start, end, regions, block, hgvs_start, hgvs_end = rows
					if regions != ivs_number:
						continue
					max_p = max(max_p, int(start), int(end))
					min_p = min(min_p, int(start), int(end))
					hgvs.regions = regions
					hgvs.blockAttr = block
				if hgvs.offset_prefix == "-":
					position = max_p - int(hgvs.offset) + 1 if strand == "+" else min_p + int(hgvs.offset) - 1
				else:
					position = min_p + int(hgvs.offset) - 1 if strand == "+" else max_p - int(hgvs.offset) + 1
				hgvs.start = position - 1
				if hgvs.mutation_type == "ins" or hgvs.mutation_type == "dup":
					hgvs.end = hgvs.start
				else:
					hgvs.end = hgvs.start + len(hgvs.ref_allele)
			elif hgvs.kind == "genomic":
				pass
			else:
				cdna_start = hgvs.cdna_start
				cdna_end = hgvs.cdna_end
				start_offset = int(hgvs.cdna_start_offset) if hgvs.cdna_start_offset else None
				end_offset = int(hgvs.cdna_end_offset) if hgvs.cdna_end_offset else None
				end_offset_prefix = hgvs.cdna_end_offset_prefix if hgvs.cdna_end_offset_prefix else None
				start_offset_prefix = hgvs.cdna_start_offset_prefix if hgvs.cdna_start_offset_prefix else None
				p_s, region_s, block_s = self._prase_cdna_pos(refseq_select, strand, cdna_start,
				                                              start_offset_prefix, start_offset) or (-1, ".", ".")
				p_e, region_e, block_e = self._prase_cdna_pos(refseq_select, strand, cdna_end,
				                                              end_offset_prefix, end_offset) or (-1, ".", ".")
				hgvs.start, hgvs.end = (int(p_s) - 1, int(p_e)) if int(p_s) < int(p_e) else (int(p_e) - 1, int(p_s))
				if hgvs.Strand == '+':
					hgvs.regions = "{0}-{1}".format(region_s, region_e) if region_s != region_e else region_s
					hgvs.blockAttr = "{0}-{1}".format(block_s, block_e) if block_s != block_e else block_s
				else:
					hgvs.regions = "{0}-{1}".format(region_e, region_s) if region_s != region_e else region_e
					hgvs.blockAttr = "{0}-{1}".format(block_e, block_s) if block_s != block_e else block_e
				if hgvs.kind == "MT":
					if hgvs.start < 0:
						hgvs.start = int(cdna_start) - 1
					if hgvs.end < 0:
						hgvs.end = int(cdna_end)
			hgvs.ref_allele, hgvs.alt_allele, hgvs.mut_statue = self.get_genomics_base(hgvs.chrom, hgvs.start,
			                                                                           hgvs.end, hgvs.ref_allele,
			                                                                           hgvs.alt_allele)
			if hgvs.mutation_type == "ins" or hgvs.mutation_type == "dup":
				hgvs.start += 1
				if hgvs.alt_allele != "." and hgvs.mutation_type == "dup":
					hgvs.ref_allele = hgvs.alt_allele
				if hgvs.ref_allele != ".":
					if hgvs.mutation_type == "dup":
						hgvs.alt_allele = hgvs.ref_allele * 2
					elif hgvs.mutation_type == "ins":
						hgvs.alt_allele = hgvs.alt_allele
						hgvs.ref_allele = hgvs.ref_allele[0]
			hgvs.end = hgvs.start + len(hgvs.ref_allele)
			if re.compile("del").search(hgvs.mutation_type):
				ref_allele_add, alt_allele_add, mut_statue_add = self.get_genomics_base(hgvs.chrom,
				                                                                        hgvs.start - 1, hgvs.start)
				hgvs.alt_allele = ref_allele_add + hgvs.alt_allele if hgvs.alt_allele != "." else ref_allele_add
				hgvs.ref_allele = ref_allele_add + hgvs.ref_allele if hgvs.ref_allele != "." else ref_allele_add
			if hgvs.mutation_type == "dup":
				ref_allele_add, alt_allele_add, mut_statue_add = self.get_genomics_base(hgvs.chrom,
				                                                                        hgvs.start - 2, hgvs.start - 1)
				hgvs.alt_allele = ref_allele_add + hgvs.ref_allele
				hgvs.ref_allele = ref_allele_add
				hgvs.start -= 1
				hgvs.end -= 1
			final = (trans_name, chgvs, hgvs.PrimaryTag, hgvs.geneSym, hgvs.protein, hgvs.Strand, hgvs.blockAttr,
			         hgvs.regions, hgvs.chrom, hgvs.start, hgvs.end, hgvs.ref_allele, hgvs.alt_allele)
			feback.add(self.hgvsname._make(final))
		return feback

	def get_genomics_base(self, chrom, start, end, refer=".", alter="."):
		chrom = str(chrom)
		if chrom.startswith("chrM"):
			chrom = 'chrM_NC_012920.1'
		if self.refer is None:
			return refer, alter, "UnKown"
		ref_base = self.refer.fetch(chrom, int(start), int(end)).upper()
		base_pairs = get_pair(ref_base)
		if not len(refer):
			refer = "."
		if not len(alter):
			alter = "."
		if "N" in refer or refer == ref_base:
			return ref_base, alter, "Match"
		elif refer == "." and alter != ".":
			return ref_base, ref_base[0] + alter, "Match"
		elif refer == "." and alter == ".":
			return ref_base, alter, "Match"
		elif base_pairs == refer:
			return ref_base, get_pair(alter), "Match"
		else:
			return ref_base, alter, "UnMatch"

	@staticmethod
	def _prase_cdna_pos(refseq, strand, pos, offset_prefix=None, offset=None):
		for rows in refseq:
			start, end, regions, block, hgvs_start, hgvs_end = rows
			utr3 = hgvs_end if strand == "-" else hgvs_start
			utr3_pos = end if strand == "-" else start
			utr5 = hgvs_start if strand == "-" else hgvs_end
			urt5_pos = start if strand == "-" else end
			if pos.startswith("*"):
				if utr3 != "*1":
					continue
				pos_r = int(pos[1:])
				if offset_prefix and offset:
					pos_r = pos_r + offset if offset_prefix == "+" else pos_r - offset
				pos = utr3_pos - (pos_r - 1) if strand == "-" else utr3_pos + (pos_r - 1)
				return pos, regions, block
			else:
				if regions.startswith("IVS") or utr3.startswith("*"):
					continue
				if int(utr3) <= int(pos) <= int(utr5):
					if strand == "+":
						pos_r = urt5_pos - (int(utr5) - int(pos))
						if offset_prefix and offset:
							pos_r = pos_r + offset if offset_prefix == "+" else pos_r - offset
					else:
						pos_r = urt5_pos + (int(utr5) - int(pos))
						if offset_prefix and offset:
							pos_r = pos_r + offset if offset_prefix == "-" else pos_r - offset
					return pos_r, regions, block

	def get_exons(self, fileout=None):
		feback = defaultdict(set)
		sql = "SELECT rv.geneSym, rv.Trans, rv.Version, ms.Regions, ms.chrom, ms.Start, ms.End FROM " \
		      "Refcontent AS ms, RefHeader AS rv WHERE ms.RefIds = rv.RefIds"
		if self.trans is not None:
			for tran in self.trans:
				condition = "  AND rv.Trans=\"%s\"" % tran
				sencur = self.dbref.execute(sql + condition).fetchall()
				if sencur is None:
					print "Sorry, Transcript %s has not been included !" % tran
					continue
				for g, t, v, r, c, s, e in sencur:
					if r.startswith("IVS"):
						continue
					tran = ".".join([t, v])
					if int(s) > int(e):
						s, e = e, s
					rec = (tran, g, r, c, str(s), str(e))
					feback[tran].add(rec)
		elif self.genes is not None:
			for gene in self.genes:
				sql += "  rv.geneSym=\"%s\"" % gene
				sencur = self.dbref.execute(sql).fetchall()
				if sencur is None:
					print "Sorry, geneSym %s has not been included !" % gene
					continue
				for g, t, v, r, c, s, e in sencur:
					if r.startswith("IVS"):
						continue
					tran = ".".join([t, v])
					if int(s) > int(e):
						s, e = e, s
					rec = (tran, g, r, c, str(s), str(e))
					feback[tran].add(rec)
		if fileout:
			fileout = os.path.abspath(fileout)
			f_out = smart_open(fileout, 'w')
			for k, v in feback.iteritems():
				f_out.writelines(["\t".join(line) + '\n' for line in v])
			f_out.close()
			return fileout
		else:
			return feback
