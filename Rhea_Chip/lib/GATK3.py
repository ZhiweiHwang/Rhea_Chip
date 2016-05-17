#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
from JavaTools import JavaTool


class GATK3(JavaTool):
	def __init__(self, jar="", javatmp=os.getcwd(), max_memory=None, **kwargs):
		if not jar:
			raise IOError("Fail to load GATK3 tools~")
		JavaTool.__init__(self, os.path.abspath(jar), tmp_path=os.path.abspath(javatmp), max_memory=max_memory)

	def RealignerTargetCreator(self, reference, in_bam, target_intervals=None, max_threads=0, knowns=None,
	                           regions=None, min_reads_cov=None, mismatch_fraction=None, **kwargs):
		in_bam = os.path.abspath(in_bam)
		prefix, _ = os.path.splitext(in_bam)
		options = ' -nt %i ' % int(max_threads) if int(max_threads) else ""
		options += " -R %s" % os.path.abspath(reference)
		options += " -I %s" % in_bam
		target_intervals = os.path.abspath(target_intervals) if target_intervals else prefix + '.realn_data.intervals'
		options += " -o %s" % target_intervals
		if knowns:
			for i in knowns.split():
				if os.path.exists(i):
					options += " --known %s" % os.path.abspath(i)
		options += " -minReads %i" % int(min_reads_cov) if min_reads_cov else ""
		options += " -mismatch %f" % float(mismatch_fraction) if mismatch_fraction else ""
		options += " -L %s" % os.path.abspath(regions) if regions else ""
		return "{0} -T RealignerTargetCreator {1}".format(self.options, options)

	def IndelRealigner(self, reference, in_bam, output=None, target_intervals=None, max_threads=0, filter_no_bases=True,
	                   knowns=None, lod_threshold=0, entropy_threshold=0, max_cons=0, max_size_for_movement=0,
	                   max_pos_move=0, max_reads_for_cons=0, max_reads_for_realignment=None, max_reads_in_memory=None,
	                   no_original_tags=0, nway_out=0, default_base_qualities=None, **kwargs):
		in_bam = os.path.abspath(in_bam)
		prefix, _ = os.path.splitext(in_bam)
		options = ' -nt %i ' % int(max_threads) if int(max_threads) else ""
		options += " -R %s" % os.path.abspath(reference)
		options += " -I %s" % in_bam
		options += " -targetIntervals %s" % os.path.abspath(target_intervals)
		output = os.path.abspath(output) if output else prefix + '.sort.realn.bam'
		options += " -o %s" % output
		options += " -filterNoBases" if filter_no_bases else ""
		if knowns:
			for i in knowns.split():
				if os.path.exists(i):
					options += " --known %s" % os.path.abspath(i)
		options += " -LOD %f" % float(lod_threshold) if lod_threshold else ""
		options += " -entropy %f" % float(entropy_threshold) if entropy_threshold else ""
		options += " --maxConsensuses %i" % int(max_cons) if max_cons else ""
		options += " -maxIsize %i" % int(max_size_for_movement) if max_size_for_movement else ""
		options += " -maxPosMove %i" % int(max_pos_move) if max_pos_move else ""
		options += " -greedy %i" % int(max_reads_for_cons) if max_reads_for_cons else ""
		options += " -maxReads %i" % int(max_reads_for_realignment) if max_reads_for_realignment else ""
		options += " -maxInMemory %i" % int(max_reads_in_memory) if max_reads_in_memory else ""
		options += " --nWayOut" if int(nway_out) else ""
		options += " -noTags" if int(no_original_tags) else ""
		options += " --defaultBaseQualities %i" % int(default_base_qualities) if default_base_qualities else ""
		return "{0} -T IndelRealigner {1}".format(self.options, options)

	def BaseRecalibrator(self, reference, in_bam, recal_data=None, knowns=None, regions=None, max_threads=0, **kwargs):
		in_bam = os.path.abspath(in_bam)
		prefix, _ = os.path.splitext(in_bam)
		options = ' -nct %i ' % int(max_threads) if int(max_threads) else ""
		options += " -R %s" % os.path.abspath(reference)
		options += " -I %s" % in_bam
		recal_data = os.path.abspath(recal_data) if recal_data else prefix + '.recal_data.grp'
		options += " -o %s" % recal_data
		if knowns:
			for i in knowns.split():
				if os.path.exists(i):
					options += " -knownSites %s" % os.path.abspath(i)
		options += " -L %s" % os.path.abspath(regions) if regions else ""
		return "{0} -T BaseRecalibrator {1}".format(self.options, options)

	def PrintReads(self, reference, in_bam, recal_data, output=None, max_threads=0, **kwargs):
		in_bam = os.path.abspath(in_bam)
		prefix, _ = os.path.splitext(in_bam)
		options = ' -nct %i ' % int(max_threads) if int(max_threads) else ""
		options += " -R %s" % os.path.abspath(reference)
		options += " -I %s" % in_bam
		options += " -BQSR %s " % os.path.abspath(recal_data)
		output = os.path.abspath(output) if output else prefix + '.recal_data.bam'
		options += " -o %s" % output
		return "{0} -T PrintReads -filterNoBases {1}".format(self.options, options)

	def HaplotypeCaller(self, reference, in_bam, output, regions=None, cov_threshold=0, max_threads=0,
	                    call_conf=30.0, emit_conf=10.0, read_filter='BadCigar', **kwargs):
		in_bam = os.path.abspath(in_bam)
		output = os.path.abspath(output)
		options = ' -nct %i ' % int(max_threads) if int(max_threads) else ""
		options += " -R %s" % os.path.abspath(reference)
		options += " -I %s" % in_bam
		options += " -L %s" % os.path.abspath(regions) if regions else " --emitRefConfidence GVCF"
		options += " -o %s" % output
		options += " -rf %s" % read_filter if read_filter else ""
		options += " -stand_call_conf %f" % float(call_conf) if float(call_conf) else ""
		options += " -stand_emit_conf %f" % float(emit_conf) if float(emit_conf) else ""
		options += " -dcov %i" % int(cov_threshold) if int(cov_threshold) else ""
		return "{0} -l INFO -T HaplotypeCaller {1} -A AlleleBalance -A HaplotypeScore".format(self.options, options)

	def UnifiedGenotyper(self, reference, in_bam, output, regions=None, cov_threshold=0, max_threads=0, call_conf=30.0,
	                     emit_conf=10.0, read_filter='BadCigar', discovery_mode="BOTH", **kwargs):
		in_bam = os.path.abspath(in_bam)
		output = os.path.abspath(output)
		options = ' -nt %i ' % int(max_threads) if int(max_threads) else ""
		options += " -R %s" % os.path.abspath(reference)
		options += " -I %s" % in_bam
		options += " -L %s" % os.path.abspath(regions) if regions else " --output_mode EMIT_ALL_SITES"
		options += " -o %s" % output
		options += " -stand_call_conf %f" % float(call_conf) if float(call_conf) else ""
		options += " -stand_emit_conf %f" % float(emit_conf) if float(emit_conf) else ""
		options += " -rf %s" % read_filter if read_filter else ""
		options += " -dcov %i" % int(cov_threshold) if int(cov_threshold) else ""
		options += " -glm %s" % discovery_mode if discovery_mode else ""
		return "{0} -l INFO -T UnifiedGenotyper {1}".format(self.options, options)

	def SelectVariants(self, reference, input_vcf, output_vcf, vartype=None, varfilter=None, max_threads=0, **kwargs):
		input_vcf = os.path.abspath(input_vcf)
		output_vcf = os.path.abspath(output_vcf)
		options = ' -nt %i' % int(max_threads) if int(max_threads) else ""
		options += " --variant %s -o %s" % (input_vcf, output_vcf)
		options += " -R %s" % os.path.abspath(reference)
		options += " -selectType \'%s\'" % vartype
		options += " -select \'%s\'" % varfilter if varfilter else ""
		return "{0} -T SelectVariants {1}".format(self.options, options)

	def VariantFiltration(self, reference, input_vcf, output_vcf, filter_expression=None,
	                      filter_tag="StandardFilter", **kwargs):
		input_vcf = os.path.abspath(input_vcf)
		output_vcf = os.path.abspath(output_vcf)
		options = " --variant %s -o %s" % (input_vcf, output_vcf)
		options += " -R %s" % os.path.abspath(reference)
		options += " --filterExpression \"%s\"" % filter_expression if filter_expression else ""
		options += " --filterName \"%s\"" % filter_tag if filter_tag else ""
		return "{0} -T VariantFiltration {1}".format(self.options, options)

	def ReadBackedPhasing(self, reference, input_bam, input_vcf, output_vcf, **kwargs):
		input_bam = os.path.abspath(input_bam)
		input_vcf = os.path.abspath(input_vcf)
		output_vcf = os.path.abspath(output_vcf)
		options = " -R %s -I %s --variant %s -o %s" % (os.path.abspath(reference), input_bam, input_vcf, output_vcf)
		return "{0} -T ReadBackedPhasing {1}".format(self.options, options)

	def CombineVariants(self, reference, input_vcf, output_vcf, max_threads=0, merge_option="UNSORTED", **kwargs):
		output_vcf = os.path.abspath(output_vcf)
		reference = os.path.abspath(reference)
		options = " -R %s -o %s" % (reference, output_vcf)
		options += " -nt %i" % int(max_threads) if int(max_threads) else ""
		options += " -genotypeMergeOptions %s" % merge_option if merge_option else ""
		for vcf in input_vcf.split():
			vcf = os.path.abspath(vcf)
			options += " --variant %s" % vcf
		return "{0} -T CombineVariants {1} ".format(self.options, options)

	def VariantRecalibrator(self, reference, input_vcf, recalfile, tranchesfile, training_set_good,
	                        training_set_bad=None, mode="SNP", **kwargs):
		input_vcf = os.path.abspath(input_vcf)
		reference = os.path.abspath(reference)
		recalfile = os.path.abspath(recalfile)
		tranchesfile = os.path.abspath(tranchesfile)
		options = " -R %s -input %s -recalFile %s -tranchesFile %s" % (reference, input_vcf, recalfile, tranchesfile)
		options += " -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an QD -mode %s " % mode
		for i in training_set_good.split():
			if os.path.exists(i):
				options += " -resource:good,known=true,training=true,truth=true %s" % os.path.abspath(i)
		if training_set_bad:
			for i in training_set_bad.split():
				if os.path.exists(i):
					options += " -resource:bad,known=true,training=true,bad=true %s" % os.path.abspath(i)
		return "{0} -T VariantRecalibrator {1} ".format(self.options, options)

	def ApplyRecalibration(self, reference, input_vcf, recalfile, tranchesfile, output_vcf,
	                       mode="SNP", ts_filter_level=99.0, **kwargs):
		input_vcf = os.path.abspath(input_vcf)
		reference = os.path.abspath(reference)
		recalfile = os.path.abspath(recalfile)
		tranchesfile = os.path.abspath(tranchesfile)
		output_vcf = os.path.abspath(output_vcf)
		options = " -R %s -input %s -recalFile %s -tranchesFile %s" % (reference, input_vcf, recalfile, tranchesfile)
		options += " -o %s  -mode %s" % (output_vcf, mode)
		options += " --ts_filter_level %f" % float(ts_filter_level) if float(ts_filter_level) else ""
		return "{0} -T ApplyRecalibration {1} ".format(self.options, options)
