#! /usr/bin/env python
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！
import json
import os
import re
from collections import namedtuple
from ConfigParser import ConfigParser

__title__ = '__init__.py'
__version__ = 'v1.0'
__author__ = 'huang'
__mtime__ = '14:27 02/22 2016'
__description__ = 'A Simple pipeline use for phrase VCF.'
__note__ = "Please specify environment variable Rhea_Chip_Home as the home of Rhea_Chip pipeline."


class AnnoVar(object):
	"""
	[general]
	database = dbsnp138
	mode = hg19
	# Reference model, option, only use in regions search, support hg19 / GRCh37
	OutputOrder = 3
	# Output Order, the smaller the value, the more output in front
	type = indexed_tsv
	# File Format, required, support indexed_vcf / indexed_tsv / sqlite3 / mysql / protalsql / tsv
	src_file = dbsnp_138.hg19.vcf.gz
	indexed_key = {"sequence" : "Chrom", "begin" : "Start", "end" : "Stop"}
	# File relative path, required in indexed_vcf / indexed_tsv / sqlite3 / tsv, right is the name used in our input
	src_config = {}
	# Sql used config, required in mysql / protalsql
	search_key = {"Chrom":{"#Chr":"string"}, "Start":{"Start":"int"}, "Stop":{"Stop":"int"}, "Refer":{"Refer":"string"}, "Call":{"Call":"string"}}
	# search key, left is the name used in our input, right is the confirm name marked in annotation file
	# left choose: Chrom / Start / Stop / Refer / Call / GeneSym / Transcript
	search_value = {"rsID":"rsID", "dbSNPAlleleFreq":"dbSNPAF"}
	# search key, left is the confirm name marked in annotation file, right is the name used in our output
	"""

	def __init__(self, config):
		self.config = os.path.abspath(config)
		handle = ConfigParser()
		handle.read(self.config)
		self.db = handle.get("general", "database")
		self.order = float(handle.get("general", "OutputOrder"))
		self.mode = handle.get("general", "mode")
		self.type = handle.get("general", "type")
		self.search_key = json.loads(handle.get("general", "search_key"))
		self.indexed_key = json.loads(handle.get("general", "indexed_key"))
		self.search_value = json.loads(handle.get("general", "search_value"))
		self.src_file = os.path.join(os.path.dirname(self.config), handle.get("general", "src_file"))
		self.src_config = json.loads(handle.get("general", "src_config"))
		self.header = list()
		self.cur = None
		self.info = namedtuple(self.db, self.search_value.values())
		if self.type in ['indexed_vcf', 'indexed_tsv']:
			import pysam
			self.handle = pysam.Tabixfile(self.src_file)
			self.header = list(self.handle.header)[-1].split("\t")
			self.fetch = self.search_from_tabix
		elif self.type == "sqlite3":
			import sqlite3
			self.handle = sqlite3.connect(self.src_file)
			self.cur = self.handle.cursor()
			self.fetch = self.search_from_sql
		elif self.type == "mysql":
			import MySQLdb
			self.handle = MySQLdb.connect(host=self.src_config['host'], user=self.src_config['user'],
			                              passwd=self.src_config['passwd'], db=self.src_config['db'],
			                              port=int(self.src_config['port']))
			self.cur = self.handle.cursor()
			self.fetch = self.search_from_sql
		elif self.type == "PostgreSQL":
			import psycopg2
			self.handle = psycopg2.connect(database=self.src_config['db'], user=self.src_config['user'],
			                               password=self.src_config['passwd'], host=self.src_config['host'],
			                               port=int(self.src_config['port']))
			self.cur = self.handle.cursor()
			self.fetch = self.search_from_sql
		elif self.type == "tsv":
			from smart_open import smart_open
			self.handle = smart_open(self.src_file)
			line = self.handle.readline()
			while not len(self.header) and line:
				rows = line.strip().split("\t")
				if self.search_key.keys()[0] in rows:
					self.header = rows
					break
				line = self.handle.readline()
			self.fetch = self.search_from_tsv

	def __del__(self):
		self.handle.close()

	@staticmethod
	def remove_base(pos, refer, alter):
		start = int(pos) - 1
		refer = refer.upper()
		alters = alter.upper()
		tmp_ref = refer
		return_set = set()
		for alter in re.split("[^ATCGN\.]+", alters):
			while len(refer) and len(alter) and refer[0] == alter[0]:
				refer = refer[1:]
				alter = alter[1:]
				start += 1
			stop = start + len(refer)
			if not len(refer) and not len(alter):
				refer = alter = tmp_ref[0]
			if not len(refer):
				refer = "."
			if not len(alter):
				alter = "."
			return_set.add((start, stop, refer, alter))
		return return_set

	def _check_result(self, chrom, start, end, **kwargs):
		nor_dict = dict()
		feback_set = set()
		if "Refer" in self.search_key and "Call" in self.search_key:
			refer_k = self.search_key["Refer"].keys()[0]
			alter_k = self.search_key["Call"].keys()[0]
		else:
			refer_k = alter_k = None
		for i, j in kwargs.iteritems():
			if i not in self.indexed_key.values() and i in self.search_key:
				for a, _ in self.search_key[i].iteritems():
					nor_dict[a] = j
		for line in self.handle.fetch(chrom, start, end):
			rows = line.strip().split("\t")
			mes = dict(zip(self.header, rows))
			check = True
			bases_level = 0
			try:
				if refer_k in mes and alter_k in mes:
					for base_set in self.remove_base(0, mes[refer_k], mes[alter_k]):
						bases_level += 1
						_, _, mes[refer_k], mes[alter_k] = base_set
						if (refer_k in nor_dict and mes[refer_k] != nor_dict[refer_k]) or \
								(alter_k in nor_dict and mes[alter_k] != nor_dict[alter_k]):
							continue
				for i in nor_dict:
					if i not in mes:
						continue
					if mes[i] != nor_dict[i]:
						raise Exception("Well, that's ok ~")
			except Exception:
				check = False
			if check is False:
				continue
			result = dict()
			for i, j in self.search_value.iteritems():
				result[j] = "."
				if i in mes:
					result[j] = mes[i]
				elif "." in i:
					first, second = i.split(".", 1)
					if first in mes:
						tmp_dict = dict()
						for info in mes[first].split(";"):
							infos = info.split("=", 1)
							tmp_dict[infos[0]] = infos[1] if len(infos) == 2 else "."
						if second in tmp_dict:
							result[j] = tmp_dict[second]
			feback_set.add(self.info(**result))
		return feback_set

	def search_from_tabix(self, **kwargs):
		feback = set()
		try:
			chrom = str(kwargs[self.indexed_key["sequence"]])
			start = int(kwargs[self.indexed_key["begin"]])
			end = int(kwargs[self.indexed_key["end"]])
			if start == end:
				start -= 1
		except KeyError:
			raise KeyError("While you use %s, you should support %s" % (self.db, self.indexed_key.values()))
		except Exception as e:
			raise Exception(e)
		if re.compile("chr", re.I).match(self.indexed_key["sequence"]):
			chrom = re.sub("^chr", "", chrom)
			if self.mode == 'hg19':
				chrom = "chr" + chrom
		if "refer" in kwargs and "alter" in kwargs:
			for base_set in self.remove_base(start + 1, kwargs["Refer"], kwargs["Call"]):
				start, end, kwargs["Refer"], kwargs["Call"] = base_set
				feback |= self._check_result(chrom, start, end, **kwargs)
		else:
			feback |= self._check_result(chrom, start, end, **kwargs)
		if not len(feback):
			result = {i: "." for i in self.search_value.values()}
			feback.add(self.info(**result))
		return feback

	def search_from_sql(self, **kwargs):
		if self.cur is None:
			raise IOError("Sorry, failed to connect to the database : %s !!" % self.db)
		feback = set()
		tmp_result = self.search_value.keys()
		sql = "SELECT %s FROM %s " % (", ".join(tmp_result), self.src_config['table'])
		lim = list()
		nor_dict = dict()
		for i, j in kwargs.iteritems():
			if i in self.search_key:
				for a, b in self.search_key[i].iteritems():
					nor_dict[a] = j
					if b == "int":
						try:
							lim.append("{0}={1}".format(a, int(j)))
						except ValueError:
							lim.append("{0}=\"{1}\"".format(a, j))
					else:
						lim.append("{0}=\"{1}\"".format(a, j))
		sql += "WHERE " + " and ".join(lim) if len(lim) else ""
		try:
			ver = self.cur.execute(sql).fetchall()
			for rows in ver:
				mes = dict(zip(tmp_result, rows))
				result = dict()
				for i, j in self.search_value.iteritems():
					result[j] = mes[i] if i in mes else "."
				feback.add(self.info(**result))
		except Exception:
			pass
		if not len(feback):
			result = dict()
			for i in self.search_value.values():
				result[i] = nor_dict[i] if i in nor_dict else "."
			feback.add(self.info(**result))
		return feback

	def search_from_tsv(self, **kwargs):
		feback = set()
		nor_dict = dict()
		for i, j in kwargs.iteritems():
			if i in self.search_key:
				for a, _ in self.search_key[i].iteritems():
					nor_dict[a] = j
		for line in self.handle.readlines():
			rows = line.strip().split("\t")
			if len(rows) != self.header:
				continue
			mes = dict(zip(self.header, rows))
			check = [mes[i] == nor_dict[i] for i in nor_dict.keys() if i in mes]
			if False in check:
				continue
			result = dict()
			for i, j in self.search_value.iteritems():
				result[j] = mes[i] if i in mes else "."
			feback.add(self.info(**result))
		if not len(feback):
			result = dict()
			for i in self.search_value.values():
				result[i] = nor_dict[i] if i in nor_dict else "."
			feback.add(self.info(**result))
		return feback
