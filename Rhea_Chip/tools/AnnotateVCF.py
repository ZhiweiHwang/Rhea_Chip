#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！
# Program       : AnnotateVCF.py
# Description   : annotate vcf's tsv result call by RheaChip .
# Dependency    : This script depend on the bgzip, tabix and samtools
#                 to be available, the whole genome fasta to be read from
#                 Rhea_Chip/db/aln_db/hg19/ named "hg19_chM_male_mask.fa", which change
#                 the original chrM of hg19 to chrM_NC_012920.1.
#                 Also the entire annotation databases to be read from Rhea_Chip/db/transdb are needed to annotate.

import logging
import os
import pysam
import re
import Rhea_Chip
import sys
import subprocess
from collections import defaultdict, namedtuple
from smart_open import smart_open
from glob import glob
from Rhea_Chip.lib.annovar import VCFParser, AnnoVar
from Rhea_Chip.tools import _chrom_valued
from Rhea_Chip.lib import snpEff
from optparse import OptionParser, OptionGroup

Rhea_Chip_Home = os.getenv('Rhea_Chip_Home') or os.path.dirname(os.path.abspath(__file__))
logger = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s: %(levelname)s: %(message)s', level=logging.INFO)
logger.info('Running: python %s' % ' '.join(sys.argv))
_TransMessages = dict()
_AnnotationDB = defaultdict(set)


class DBAnno(object):
	def __init__(self, annotationdb=_AnnotationDB):
		self._dbhandles = list()
		dbtitle = defaultdict(float)
		for dbconfigset in annotationdb.values():
			for dbconfig in dbconfigset:
				dbs = AnnoVar(dbconfig)
				t = "\t".join([str(i) for i in sorted(dbs.search_value.values())])
				dbtitle[t] = dbs.order
				self._dbhandles.append(dbs)
		self.dbtitle = "\t".join([j for j in sorted(dbtitle.keys(), key=lambda x: dbtitle[x])])

	def __del__(self):
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
			except KeyError as err:
				logger.info(err)
				continue
			except Exception:
				continue
		return dbinfo


class VCFReader(object):
	def __init__(self, vcf, reference):
		self.vcf = VCFParser.Reader(vcf)
		try:
			pattern = re.compile("VCFv(\d+)\.(\d+)").match(self.vcf.metadata['fileformat']).groups()
			version = (int(pattern[0]), int(pattern[1]))
			if version < (4, 0):
				raise AttributeError("Only handle VCF 4.0 or later versions !")
			if version != (4, 1):
				logger.info("As your fileformat version is %s, assuming it as VCFv4.1" % ".".join(pattern))
		except Exception:
			raise AttributeError("Sorry, %s is not a regular VCF file" % vcf)
		self.refer = pysam.FastaFile(reference)
		self.DBAnno = DBAnno(_AnnotationDB)

	def __del__(self):
		self.vcf.__del__()
		self.refer.close()

	def get_snpeff_info(self, snpeff_info):
		global _TransMessages
		try:
			eff_info = re.split("\s*\|\s*", re.compile("Functional annotations:\s+'([^']+)'\s+").
			                    match(self.vcf.infos['ANN'].desc).group(1))
		except Exception:
			raise AttributeError("Sorry, we can't get anno message from %s" % self.vcf.filename)
		eff_info_dict = defaultdict(set)
		for eff in snpeff_info.split(","):
			rows = re.split("\s*\|\s*", eff)
			eff_dict = defaultdict(str, zip(eff_info, rows))
			alter = eff_dict['Allele'] or "."
			if alter == ".":
				continue
			functions = eff_dict['Annotation'] or "."
			genesym = eff_dict['Gene_Name'] or "Intergenic"
			chgvs = eff_dict['HGVS.c'] or "."
			phgvs = eff_dict['HGVS.p'] or "."
			trans = eff_dict['Feature_ID'].split(".")[0] or "."
			if eff_dict['Feature_Type'] != 'transcript' and trans not in _TransMessages:
				continue
			if len(_TransMessages) and trans not in _TransMessages and trans.startswith("NM_"):
				continue
			trans_bio = eff_dict['Transcript_BioType'] or "."
			impact = eff_dict['Annotation_Impact'] or "."
			try:
				exon_num = re.compile("(\d+)/\d+").match(eff_dict['Rank']).group(1) or "."
			except Exception:
				exon_num = ""
			if re.compile(".\.-").match(chgvs) and not exon_num:
				exon_num = "PROM"
			elif re.compile(".\.\d+[+-]\d+").match(chgvs):
				exon_num = "IVS" + exon_num if exon_num else "intergenic"
			else:
				exon_num = "EX" + exon_num if exon_num else "intergenic"
			cdna_pos = re.compile("(\d+)/\d+").match(eff_dict['cDNA.pos / cDNA.length'])
			cdna_pos = cdna_pos.group(1) if cdna_pos else "."
			aa_pos = re.compile("(\d+)/\d+").match(eff_dict['AA.pos / AA.length'])
			aa_pos = aa_pos.group(1) if aa_pos else "."
			cds_pos = re.compile("(\d+)/\d+").match(eff_dict['CDS.pos / CDS.length'])
			cds_pos = cds_pos.group(1) if cds_pos else "."
			anno_tag = eff_dict['ERRORS / WARNINGS / INFO'] or "."
			if trans in _TransMessages:
				gene_id, trans, protein_id, stand, primary = _TransMessages[trans]
			else:
				gene_id = eff_dict['Gene_ID'] or "."
				protein_id = "."
				stand = "."
				primary = "Y"
			eff_info_dict[alter].add((genesym, gene_id, trans, trans_bio, chgvs, protein_id, phgvs, stand, primary,
			                          functions, impact, exon_num, cdna_pos, aa_pos, cds_pos, anno_tag))
		return eff_info_dict

	def parse(self, outdir=os.getcwd(), kickout_function=None, follow_function=None):
		logger.info("General VCF Phrasing begin !!!")
		self.vcf.parse_metainfo()
		titles = ["Chrom", "Start", "Stop", "Refer", "Call", "Zygosity", "VarType", "Filter", "ADepth", "ARatio",
		          "PL", "NeighborGID", "PhasedGID", "MutationName", "GeneSym", "EntrezGeneID", "Transcript",
		          "TransBioType", "cHGVS", "Protein", "pHGVS", "Strand", "PrimaryTag", "FunctionName", "Impact",
		          "ExInID", "cDNAPos", "AAPos", "CDSPos", "AnnoTag", "StandardMutation"]
		annotation = namedtuple("Annotation", titles)
		dbtitle = self.DBAnno.dbtitle.split("\t")
		titles.extend(dbtitle)
		kickouts = kickout_function.split(",") if kickout_function else list()
		follows = follow_function.split(",") if follow_function else None
		for name, samples in self.vcf.readlines.iteritems():
			logger.info("Start to annotate sample %s" % name)
			message_num = 0
			result_out = list()
			for messages in samples:
				chrom = messages.Chrom
				pos = int(messages.Pos)
				gt = re.compile("(\d)[|/](\d)").match(messages.GT)
				if not gt:
					continue
				message_num += 1
				if message_num % 250 == 0:
					logger.info("Complete " + str(message_num) + " articles")
				l1, l2 = int(gt.group(1)), int(gt.group(2))
				refer = messages.Ref
				alter_1 = messages.Alter[l1 - 1] if l1 > 0 else refer
				alter_2 = messages.Alter[l2 - 1] if l2 > 0 else refer
				alter_out_set = set()
				eff_info_dict = self.get_snpeff_info(messages.Info['ANN'])
				try:
					nb_id = messages.NB
				except AttributeError:
					nb_id = "."
				try:
					pb_id = messages.PB
				except AttributeError:
					pb_id = "."
				filter_tag = messages.Filter
				if l1 == l2:
					zygosity = "hom-ref" if l1 == 0 else "hom-alt"
					alter_out_set.add(alter_1)
				else:
					zygosity = "het-alt"
					if l1 != 0:
						alter_out_set.add(alter_1)
					if l2 != 0:
						alter_out_set.add(alter_2)
				for alters in alter_out_set:
					offset_s, offset_e, nor_refer, nor_alter, vartype, close_anno = \
						self.closet_anno(chrom, pos, refer, alters)
					start = pos + offset_s
					stop = pos + + offset_e
					for info_rows in eff_info_dict[alters]:
						if not len(info_rows):
							info_rows = (".",) * 16
						genesym, gene_id, trans, trans_bio, chgvs, protein_id, phgvs, strand, primary, \
						functions, impact, exon_num, cdna_pos, aa_pos, cds_pos, anno_tag = info_rows
						functions = re.sub("_variant", "", functions)
						if (functions in kickouts) or (follows and functions not in follow_function):
							continue
						if chgvs != ".":
							mutation_name = "{0}({1}): {2}".format(trans, genesym, chgvs)
							mutation_name += " (%s)" % phgvs if phgvs != "." else ""
						else:
							mutation_name = "."
						try:
							alle_depth = int(messages.AD.split(",")[l2])
							try:
								base_depth = int(messages.DP)
								alle_radio = round(float(alle_depth) / base_depth, 2)
							except Exception:
								alle_radio = "."
						except Exception:
							alle_depth = "."
							alle_radio = "."
						try:
							phred = ",".join(messages.PL.split(",")[3 * (l2 - 1): 3 * l2])
						except Exception:
							phred = "."
						varinfo = [chrom, start, stop, nor_refer, nor_alter, zygosity, vartype, filter_tag, alle_depth,
						           alle_radio, phred, nb_id, pb_id, mutation_name, genesym, gene_id, trans, trans_bio,
						           chgvs, protein_id, phgvs, strand, primary, functions, impact, exon_num, cdna_pos,
						           aa_pos, cds_pos, anno_tag, close_anno]
						variation = annotation._make(varinfo)
						dbinfo = self.DBAnno.dbanno(variation)
						for i in dbtitle:
							if i in dbinfo:
								varinfo.append("|".join(dbinfo[i]))
							else:
								varinfo.append(".")
						result_out.append(varinfo)
			sample_out = os.path.join(outdir, "%s.anno.tsv" % name)
			f_out = open(sample_out, "w")
			f_out.write("#" + "\t".join(titles) + '\n')
			for anno_message in sorted(result_out, key=lambda x: (_chrom_valued(x[0]), int(x[1]), int(x[2]))):
				f_out.write("\t".join(map(str, anno_message)) + '\n')
			f_out.close()
			fileout = pysam.tabix_index(sample_out, seq_col=0, start_col=1, end_col=2, force=True)
			logger.info(
				"Sample {0} [{1}]: Annotation completed, total {2} articles !".format(name, fileout, message_num))

	@staticmethod
	def get_repeat(mainstring):
		strarr = str()
		count = int()
		mainstring = mainstring.upper()
		for i in range(0, len(mainstring)):
			strarr += mainstring[i]
			splitlist = mainstring.split(strarr)
			count = 0
			for j in splitlist:
				if j == '':
					count += 1
			if count == len(splitlist):
				break
		return strarr, count - 1

	def check_repeat(self, chrom, pos, bases):
		start = int(pos)
		bases, _ = self.get_repeat(bases)
		stop = start + len(bases)
		begin = max(start - 100, len(bases))
		end = min(stop + 100, self.refer.get_reference_length(chrom))
		refers = self.refer.fetch(chrom, begin, end).upper()
		offset_s = start - begin
		offset_e = stop - begin - len(bases)
		while offset_s > 0 and refers[offset_s - len(bases): offset_s] == bases:
			offset_s -= len(bases)
		while offset_e < len(refers) and refers[offset_e: offset_e + len(bases)] == bases:
			offset_e += len(bases)
		return begin + offset_s, (offset_e - offset_s) / len(bases), bases

	def closet_anno(self, chrom, start, refer, alter, base_level=1):
		offset = 0
		tmp_ref = refer
		start = int(start) - base_level
		while len(refer) and len(alter) and refer[0] == alter[0]:
			refer = refer[1:]
			alter = alter[1:]
			offset += 1
		start += offset
		if refer == alter:
			vartype = "ref"
			refer = alter = tmp_ref[0]
			offset_s = offset - 2
			r_start, repeat, bases = start, 0, refer
			result = "%s=" % refer
		elif len(refer) == len(alter) == 1:
			vartype = "snv"
			offset_s = offset - 1
			r_start, repeat, bases = start, 0, refer
			result = "%s>%s" % (refer, alter)
		elif len(refer) == 0:
			vartype = "ins"
			refer = "."
			offset_s = offset - 1
			r_start, repeat, bases = self.check_repeat(chrom, start, alter)
			if repeat > 0:
				aft = repeat + len(alter) / len(bases)
				result = "dup%s" % bases if float(aft) / repeat == 2.0 else "%s(%d>%d)" % (bases, repeat, int(aft))
			else:
				result = "ins%s" % alter
		elif len(alter) == 0:
			vartype = "del"
			alter = "."
			offset_s = offset - 1
			r_start, repeat, bases = self.check_repeat(chrom, start, refer)
			if repeat > 0:
				aft = repeat - float(len(refer)) / len(bases)
				result = "%s(%d>%d)" % (bases, repeat, aft) if aft > 0 else "del%s" % bases
			else:
				result = "del%s" % refer
		else:
			vartype = "delins"
			offset_s = offset - 1
			r_start, repeat, bases = start, 0, refer
			result = "del%sins%s" % (refer, alter)
		offset_e = offset_s + len(refer) if vartype != "ins" else offset_s
		if len(refer) > 50 or len(alter) > 50:
			vartype = "sv"
		chrom = re.sub("^chr", "", chrom)
		try:
			refer_chr = "NC_{:0>6}".format(int(chrom))
		except ValueError:
			if re.compile("m", re.I).match(chrom):
				refer_chr = 'NC_012920.1'
			elif re.compile("x", re.I).match(chrom):
				refer_chr = "NC_000023"
			elif re.compile("y", re.I).match(chrom):
				refer_chr = "NC_000024"
			else:
				refer_chr = 'chr' + chrom
		std_out = "." if repeat == 0 else "{0}:g.{1}{2}".format(refer_chr, r_start, result)
		return offset_s, offset_e, refer, alter, vartype, std_out


def main():
	global _TransMessages, _AnnotationDB
	usage = 'Usage: %prog [-h] [--version] --vcf [VCF file] --config [panel.conf] [options]'
	description = Rhea_Chip.__description__
	author = Rhea_Chip.__author__
	version = Rhea_Chip.__version__
	parser = OptionParser(usage=usage, version=version, description=description, epilog=author)
	expect = OptionGroup(parser, 'Expected arguments', 'Caution: These parameters are necessary for depthQC.')
	expect.add_option('-v', '--vcf', metavar='FILE', dest='vcf', help='vcf file')
	parser.add_option_group(expect)
	optinal = OptionGroup(parser, 'Optional arguments', 'Caution: If you do not set these parameters in addition,'
	                                                    ' AnnotateVCF will select the default value.')
	optinal.add_option('-o', '--outdir', dest='outdir', default=os.getcwd(),
	                   help='The output dir [ default: %s ]' % os.getcwd())
	optinal.add_option('-d', '--database', dest='dbdir', default=os.path.join(Rhea_Chip_Home, "db"),
	                   help='The database dir [ default: %s ]' % os.path.join(Rhea_Chip_Home, "db"))
	optinal.add_option('-r', '--reference', dest='reference',
	                   default=os.path.join(Rhea_Chip_Home, "db", "aln_db", "hg19", "hg19_chM_male_mask.fa"),
	                   help='The human reference')
	optinal.add_option('--usedb', dest='usedb', default=None,
	                   help='Definition annotation database required. separated by commas, default: all')
	optinal.add_option('--ignore', dest='ignoredb', default=None,
	                   help='Define the database you want to ignore, separated by commas, default: None')
	optinal.add_option('-t', '--trans', dest='transfile', default=None,
	                   help='Only use the transcripts in this file. Format: One transcript ID per line.')
	optinal.add_option('-g', '--gene', dest='genesfile', default=None,
	                   help='Only use the geneSym in this file. Format: One gene name per line.')
	optinal.add_option('--kickout_function', dest='kickout_function', default="no-change,.",
	                   help='kickout function, eparated by commas')
	optinal.add_option('--follow_function', dest='follow_function', default=None,
	                   help='Special focus function, eparated by commas')
	parser.add_option_group(optinal)
	(options, args) = parser.parse_args()
	if not options.vcf:
		parser.print_help()
		return 'Expected arguments lost !!!'
	vcf = os.path.abspath(options.vcf)
	outdir = os.path.abspath(options.outdir)
	databases = os.path.abspath(options.dbdir)
	usedb = options.usedb.split(",") if options.usedb else list()
	ignoredb = options.ignoredb.split(",") if options.ignoredb else list()
	transfile = os.path.abspath(options.transfile) if options.transfile else None
	genesfile = os.path.abspath(options.genesfile) if options.genesfile else None
	kickout_function = options.kickout_function
	follow_function = options.follow_function
	reference = os.path.abspath(options.reference)
	for db in glob(os.path.join(databases, "*", "*.db.config")):
		db = os.path.abspath(db)
		dbname = os.path.basename(os.path.dirname(db))
		if (dbname in ignoredb) or (len(usedb) and dbname not in usedb):
			continue
		_AnnotationDB[dbname].add(db)
	vcf_r = VCFReader(vcf, reference)
	if "ANN" not in vcf_r.vcf.infos:
		logger.info("Before we get started, we need to take some processing on the VCF !!!")
		eff = os.path.join(outdir, "snpEffAnno." + os.path.basename(vcf))
		snpeff = snpEff(os.path.join(Rhea_Chip_Home, "tools", "snpEff", "snpEff.jar"), javatmp=outdir)
		command = snpeff.VariantAnnotator(vcf, output=eff)
		_ = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		vcf_r = VCFReader(eff, reference)
	transdb = os.path.join(databases, "transdb", "transdb.config")
	transdb_handle = AnnoVar(transdb)
	transset = set()
	logger.info("Preparing transcript database !!!")
	if transfile:
		with smart_open(transfile) as f_in:
			for line in f_in:
				trans = line.strip().split()[0]
				names = trans.split(".")
				if len(names) < 1:
					continue
				transset |= transdb_handle.fetch(Trans=names[0]) if len(names) == 1 else \
					transdb_handle.fetch(Trans=names[0], TransVersion=names[1])
	elif genesfile:
		with smart_open(transfile) as f_in:
			for line in f_in:
				names = line.strip().split()[0]
				if len(names) < 1:
					continue
				transset |= transdb_handle.fetch(GeneSym=names[0])
	else:
		transset = transdb_handle.fetch(primary="Y")
	for hgvs in transset:
		trans = ".".join([hgvs.Trans, hgvs.TransVersion]) if hgvs.TransVersion != "." else hgvs.Trans
		_TransMessages[hgvs.Trans] = [hgvs.GeneID, trans, hgvs.protein_id, hgvs.Strand, hgvs.primaryTag]
	vcf_r.parse(outdir=outdir, kickout_function=kickout_function, follow_function=follow_function)


if __name__ == '__main__':
	sys.exit(main())
