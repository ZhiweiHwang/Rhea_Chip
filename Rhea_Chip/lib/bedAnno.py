#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
import sys
from glob import glob
from collections import defaultdict
from Rhea_Chip.lib.hgvs import HGVS
from Rhea_Chip.lib.annovar import AnnoVar
from smart_open import smart_open

Rhea_Chip_Home = os.getenv('Rhea_Chip_Home') or os.path.dirname(os.path.abspath(__file__))
database_dir = os.path.join(Rhea_Chip_Home, "db")

class BedAnno(object):
	def __init__(self, reference, trans=None, genes=None):
		self.refdb = os.path.join(database_dir, "transdb", "ncbi_anno_rel104.dbref.db")
		reference = os.path.abspath(reference)
		self.HGVS = HGVS(self.refdb, reference, trans, genes)
		self._dbhandles = list()
		dbtitle = defaultdict(float)
		for f in glob(os.path.join(database_dir, "*", "*.bedanno.config")):
			dbs = AnnoVar(os.path.abspath(f))
			t = "\t".join([str(i) for i in sorted(dbs.search_value.values())])
			dbtitle[t] = dbs.order
			self._dbhandles.append(dbs)
		self.dbtitle = "\t".join([j for j in sorted(dbtitle.keys(), key=lambda x: dbtitle[x])])

	def __del__(self):
		self.HGVS.__del__()
		for dbs in self._dbhandles:
			dbs.__del__()

	def dbanno(self, variation):
		dbinfo = defaultdict(set)
		for dbs in self._dbhandles:
			try:
				infos = dbs.fetch(**variation.__dict__)
				for info in infos:
					for k, v in info.__dict__.iteritems():
						if v != ".":
							dbinfo[k].add(v)
			except Exception:
				continue
		return dbinfo

	def bedanno(self, bedfile, fileout=sys.stdout):
		if not (bedfile and os.path.isfile(bedfile)):
			raise IOError("Fail to load file: %s" % bedfile)
		regions = smart_open(bedfile).readlines()
		f_out = smart_open(fileout, 'w')
		titles = "Transcript\tgeneSym\tgeneSym\tcHGVS\tProtein\tStand\tExonRegions"
		for line in regions:
			rows = line.strip().split("\t")
			if len(rows) < 3:
				f_out.write(line)
			try:
				chrom = str(rows[0])
				start = int(rows[1])
				stop = int(rows[2])
			except ValueError:
				if len(self.dbtitle):
					f_out.write("\t".join(rows) + '\t'.join([titles, self.dbtitle]) + '\n')
				else:
					f_out.write("\t".join(rows) + '\t' + titles + '\n')
				continue
			anno_m = self.HGVS.annobed(chrom, start, stop)
			for anno in anno_m:
				dbinfo = self.dbanno(anno)
				rows[0] = str(anno.Chrom)
				rows[1] = str(anno.Start)
				rows[2] = str(anno.Stop)
				trans = str(anno.Transcript)
				gene = str(anno.geneSym)
				chgvs = str(anno.cHgvs)
				protein = str(anno.Protein)
				stand = str(anno.Stand)
				exons = str(anno.ExonRegions)
				anno_message = [str(i) for i in rows + [trans, gene, chgvs, protein, stand, exons]]
				for i in self.dbtitle.split("\t"):
					if i in dbinfo:
						anno_message.append(";".join(dbinfo[i]))
					else:
						anno_message.append(".")
				f_out.write("\t".join(anno_message) + '\n')
		f_out.close()