#! /usr/bin/env python
# -*- coding: utf-8 -*-
# This file is copied from https://github.com/jamescasbon/PyVCF

import os
import re
import pysam
from smart_open import smart_open
from itertools import izip
from collections import OrderedDict, namedtuple, defaultdict

# Spec is a bit weak on which metadata lines are singular, like fileformat
# and which can have repeats, like contig
SINGULAR_METADATA = ['fileformat', 'fileDate', 'reference']

# Conversion between value in file and Python value
field_counts = {
	'.': None,  # Unknown number of values
	'A': -1,  # Equal to the number of alternate alleles in a given record
	'G': -2,  # Equal to the number of genotypes in a given record
	'R': -3,  # Equal to the number of alleles including reference in a given record
}

_Info = namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
_Filter = namedtuple('Filter', ['id', 'desc'])
_Alt = namedtuple('Alt', ['id', 'desc'])
_Format = namedtuple('Format', ['id', 'num', 'type', 'desc'])
_SampleInfo = namedtuple('SampleInfo', ['samples', 'gt_bases', 'gt_types', 'gt_phases'])
_Contig = namedtuple('Contig', ['id', 'length'])


class VcfMetadataParser(object):
	"""Parse the metadat in the header of a VCF file."""

	def __init__(self):
		super(VcfMetadataParser, self).__init__()
		self.info_pattern = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),\s*
            Number=(?P<number>-?\d+|\.|[AGR]),\s*
            Type=(?P<type>Integer|Float|Flag|Character|String),\s*
            Description="(?P<desc>[^"]*)"
            (?:,\s*Source="(?P<source>[^"]*)")?
            (?:,\s*Version="?(?P<version>[^"]*)"?)?
            >''', re.VERBOSE)
		self.filter_pattern = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),\s*
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
		self.alt_pattern = re.compile(r'''\#\#ALT=<
            ID=(?P<id>[^,]+),\s*
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
		self.format_pattern = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),\s*
            Number=(?P<number>-?\d+|\.|[AGR]),\s*
            Type=(?P<type>.+),\s*
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)
		self.contig_pattern = re.compile(r'''\#\#contig=<
            ID=(?P<id>[^>,]+)
            (,.*length=(?P<length>-?\d+))?
            .*
            >''', re.VERBOSE)
		self.meta_pattern = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')

	@staticmethod
	def vcf_field_count(num_str):
		"""Cast vcf header numbers to integer or None
		:param num_str: number or string
		"""
		if num_str is None:
			return None
		elif num_str not in field_counts:
			# Fixed, specified number
			return int(num_str)
		else:
			return field_counts[num_str]

	def read_info(self, info_string):
		"""Read a meta-information INFO line.
		:param info_string: VCF nfo
		"""
		match = self.info_pattern.match(info_string)
		if not match:
			raise SyntaxError(
				"One of the INFO lines is malformed: %s" % info_string)

		num = self.vcf_field_count(match.group('number'))

		info = _Info(match.group('id'), num,
		             match.group('type'), match.group('desc'),
		             match.group('source'), match.group('version'))

		return match.group('id'), info

	def read_filter(self, filter_string):
		"""Read a meta-information FILTER line.
		:param filter_string: filter string
		"""
		match = self.filter_pattern.match(filter_string)
		if not match:
			raise SyntaxError(
				"One of the FILTER lines is malformed: %s" % filter_string)

		filt = _Filter(match.group('id'), match.group('desc'))

		return match.group('id'), filt

	def read_alt(self, alt_string):
		"""Read a meta-information ALTline.
		:param alt_string: alter string
		"""
		match = self.alt_pattern.match(alt_string)
		if not match:
			raise SyntaxError(
				"One of the FILTER lines is malformed: %s" % alt_string)

		alt = _Alt(match.group('id'), match.group('desc'))

		return match.group('id'), alt

	def read_format(self, format_string):
		"""Read a meta-information FORMAT line.
		:param format_string: format string
		"""
		match = self.format_pattern.match(format_string)
		if not match:
			raise SyntaxError(
				"One of the FORMAT lines is malformed: %s" % format_string)

		num = self.vcf_field_count(match.group('number'))

		form = _Format(match.group('id'), num,
		               match.group('type'), match.group('desc'))

		return match.group('id'), form

	def read_contig(self, contig_string):
		"""Read a meta-contigrmation INFO line.
		:param contig_string: chrom string
		"""
		match = self.contig_pattern.match(contig_string)
		if not match:
			raise SyntaxError(
				"One of the contig lines is malformed: %s" % contig_string)
		length = self.vcf_field_count(match.group('length'))
		contig = _Contig(match.group('id'), length)
		return match.group('id'), contig

	@staticmethod
	def read_meta_hash(meta_string):
		# assert re.match("##.+=<", meta_string)
		"""
		:param meta_string: meta string
		:return: key value
		"""
		items = meta_string.split('=', 1)
		# Removing initial hash marks
		key = items[0].lstrip('#')
		# N.B., items can have quoted values, so cannot just split on comma
		val = OrderedDict()
		state = 0
		k = ''
		v = ''
		for c in items[1].strip('[<>]'):

			if state == 0:  # reading item key
				if c == '=':
					state = 1  # end of key, start reading value
				else:
					k += c  # extend key
			elif state == 1:  # reading item value
				if v == '' and c == '"':
					v += c  # include quote mark in value
					state = 2  # start reading quoted value
				elif c == ',':
					val[k] = v  # store parsed item
					state = 0  # read next key
					k = ''
					v = ''
				else:
					v += c
			elif state == 2:  # reading quoted item value
				if c == '"':
					v += c  # include quote mark in value
					state = 1  # end quoting
				else:
					v += c
		if k != '':
			val[k] = v
		return key, val

	def read_meta(self, meta_string):
		if re.match("##.+=<", meta_string):
			return self.read_meta_hash(meta_string)
		match = self.meta_pattern.match(meta_string)
		if not match:
			# Spec only allows key=value, but we try to be liberal and
			# interpret anything else as key=none (and all values are parsed
			# as strings).
			return meta_string.lstrip('#'), 'none'
		return match.group('key'), match.group('val')


class Reader(object):
	""" Reader for a VCF v 4.0 file, an iterator returning ``_Record objects`` """

	def __init__(self, vcf, chrom=None, start=None, end=None):
		""" Create a new Reader for a VCF file.

			You must specify either fsock (stream) or filename.  Gzipped streams
			or files are attempted to be recogized by the file extension, or gzipped
			can be forced with ``compressed=True``

			'prepend_chr=True' will put 'chr' before all the CHROM values, useful
			for different sources.

			'strict_whitespace=True' will split records on tabs only (as with VCF
			spec) which allows you to parse files with spaces in the sample names.
		"""
		super(Reader, self).__init__()
		self.filename = os.path.basename(vcf)
		self._row_pattern = re.compile("\t")
		self._alt_pattern = re.compile('[\[\]]')
		self._tabix = pysam.Tabixfile(vcf)
		self._reader = smart_open(vcf, 'rb')
		self.reader = (line.strip() for line in self._reader if line.strip()) if not chrom \
			else list(self.fetch(chrom, start, end))
		self.metadata = None
		self.infos = None
		self.filters = None
		self.alts = None
		self.formats = None
		self.contigs = None
		self.samples = None
		self._sample_indexes = None
		self._header_lines = list()
		self._column_headers = list()
		self.readlines = defaultdict(list)
		self._parse_header()

	def __iter__(self):
		return self.readlines

	def __del__(self):
		self._tabix.close()
		self._reader.close()

	def _parse_header(self):
		"""Parse the information stored in the metainfo of the VCF.

		The end user shouldn't have to use this.  She can access the metainfo
		directly with ``self.metadata``."""
		for attr in ('metadata', 'infos', 'filters', 'contigs', 'alts', 'formats'):
			setattr(self, attr, OrderedDict())
		parser = VcfMetadataParser()
		line = next(self.reader)
		while line.startswith('##'):
			self._header_lines.append(line)
			if line.startswith('##INFO'):
				key, val = parser.read_info(line)
				self.infos[key] = val
			elif line.startswith('##FILTER'):
				key, val = parser.read_filter(line)
				self.filters[key] = val
			elif line.startswith('##ALT'):
				key, val = parser.read_alt(line)
				self.alts[key] = val
			elif line.startswith('##FORMAT'):
				key, val = parser.read_format(line)
				self.formats[key] = val
			elif line.startswith('##contig'):
				key, val = parser.read_contig(line)
				self.contigs[key] = val
			else:
				key, val = parser.read_meta(line)
				if key in SINGULAR_METADATA:
					self.metadata[key] = val
				else:
					if key not in self.metadata:
						self.metadata[key] = []
					self.metadata[key].append(val)
			line = next(self.reader)

		fields = self._row_pattern.split(line[1:])
		self._column_headers = namedtuple("Sample", fields[:9])
		self.samples = fields[9:]
		self._sample_indexes = dict([(x, i) for (i, x) in enumerate(self.samples)])

	def fetch(self, chrom, start=None, end=None):
		if not self.filename:
			raise Exception('Please provide a filename (or a "normal" fsock)')
		reader = self._tabix.fetch(chrom, start, end)
		return reader

	def parse_metainfo(self, filter_ref=True):
		for line in self.reader:
			if line.startswith("#"):
				continue
			rows = line.split("\t")
			chrom = "chr" + re.sub("^chr", "", rows[0])
			pos = int(rows[1])
			ids = rows[2]
			ref = rows[3]
			alters = rows[4].split(",")
			try:
				qual = int(rows[5])
			except ValueError:
				try:
					qual = float(rows[5])
				except ValueError:
					qual = "."
			filt = rows[6]
			info = self._parse_info(rows[7])
			fmt = dict((y, x) for (x, y) in enumerate(rows[8].split(":")))
			for name, sample in izip(self.samples, rows[9:]):
				sampdat = ["."] * len(fmt)
				samples = sample.split(":")
				for i, j in enumerate(samples):
					sampdat[i] = j
				if filter_ref:
					if sampdat[fmt['GT']] == "." or re.compile("0[/|]0").match(sampdat[fmt['GT']]):
						continue
				message = namedtuple("Message", ["Chrom", "Pos", "ID", "Ref", "Alter",
				                                 "Qual", "Filter", "Info"] + sorted(fmt.keys(), key=lambda s: fmt[s]))
				messages = [chrom, pos, ids, ref, alters, qual, filt, info] + sampdat
				self.readlines[name].append(message._make(messages))

	@staticmethod
	def _parse_info(info_str):
		retdict = OrderedDict()
		if info_str == '.':
			return retdict
		entries = info_str.split(';')
		for entry in entries:
			entry = entry.split('=', 1)
			retdict[entry[0]] = entry[1] if len(entry) == 2 else "."
		return retdict
