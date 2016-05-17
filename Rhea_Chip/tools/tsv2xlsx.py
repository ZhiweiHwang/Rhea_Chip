#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
import sys
import re
import xlsxwriter
from collections import defaultdict
from smart_open import smart_open
from optparse import OptionParser


class Tsv2XlsX(object):
	def __init__(self, excel_out, f_format=None):
		self.excel = os.path.abspath(re.sub("\.xls.?$", "", excel_out, flags=re.I)) + '.xlsx'
		self._handle = xlsxwriter.Workbook(self.excel)
		self.OddRowFormat = self._handle.add_format(
			{"font": 'Arial', "size": 11, "bold": 1, "align": 'left', "valign": 'top', "text_wrap": 0})
		self.OddRowFormat.set_border(1)
		self.OddRowFormat.set_bg_color('#FFCC99')
		self.EvenRowFormat = self._handle.add_format(
			{"font": 'Arial', "size": 11, "bold": 1, "align": 'left', "valign": 'top', "text_wrap": 0})
		self.EvenRowFormat.set_border(1)
		self.EvenRowFormat.set_bg_color('#CCFFFF')
		self.sheet_format = self.parse_format_file(f_format) if f_format else dict()

	def __del__(self):
		self._handle.close()

	@staticmethod
	def parse_format_file(format_file=None):
		formact_dict = defaultdict(dict)
		if format_file is None or not os.path.isfile(format_file):
			return formact_dict
		try:
			varsheet_formact = smart_open(format_file)
			varsheet_lines = varsheet_formact.readlines()
			varsheet_title = varsheet_lines.pop(0).strip() if varsheet_lines[0].startswith(
				"#title") else "#title\twidth\thidden\tlevel"
			vs_title = {j: i for i, j in enumerate(varsheet_title.split("\t"))}
			for varline in varsheet_lines:
				vs_rows = varline.strip().split("\t")
				if len(vs_rows) != len(vs_title):
					continue
				k_w = re.sub("^\s+|\s+$", "", vs_rows[vs_title["#title"]].lower())
				formact_dict[k_w] = {i: int(vs_rows[j]) for i, j in vs_title.iteritems() if i != "#title"}
			varsheet_formact.close()
			return formact_dict
		except Exception:
			return formact_dict

	def tsv2xlsx(self, f_in, sheet=1):
		if not os.path.isfile(f_in):
			raise IOError("Fail to load file: %s" % f_in)
		path, name = os.path.split(f_in)
		name = ".".join(name.split('.')[1:-1])
		if len(name) > 20 or len(name) < 1:
			name = 'sheet' + str(sheet)
		worksheet = self._handle.add_worksheet(name)
		row = 0
		files = smart_open(f_in).readlines()
		titles = {str(j).lower(): int(i) for i, j in enumerate(files[0].strip().split("\t"))}
		for line in files:
			rows = line.strip().split('\t')
			row += 1
			row_name = 'A' + str(row)
			if row % 2:
				worksheet.write_row(row_name, rows, self.EvenRowFormat)
			else:
				worksheet.write_row(row_name, rows, self.OddRowFormat)
		for colname, colformat in self.sheet_format.iteritems():
			if colname not in titles:
				continue
			colnumber = int(titles[colname])
			weight = colformat["width"] or None
			otherdict = {n: m or None for (n, m) in colformat.iteritems() if n != "width"}
			worksheet.set_column(colnumber, colnumber, weight, None, otherdict)


def main():
	usage = "python %s -o <excel_out> -r <tsv_format1> <tsv_in1> <tsv_in2> ..." % sys.argv[0]
	parser = OptionParser(usage=usage)
	parser.add_option('-o', '--out', metavar='FILE', dest='output', help='The output excel file')
	parser.add_option('-r', '--format', metavar='FILE', dest='formats', help='The output format file', default=None)
	options, args = parser.parse_args()
	if not options.output:
		parser.print_help()
		return 'Expected arguments lost !!!'
	tsv_format = os.path.abspath(options.formats) if options.formats else None
	xlsx_out = Tsv2XlsX(options.output, tsv_format)
	sheet = 1
	for f_in in args:
		tsv = os.path.abspath(f_in)
		try:
			xlsx_out.tsv2xlsx(tsv, sheet)
		except IOError:
			continue
		sheet += 1
	xlsx_out.__del__()


if __name__ == '__main__':
	sys.exit(main())
