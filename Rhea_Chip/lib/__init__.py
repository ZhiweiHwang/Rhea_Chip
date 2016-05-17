#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

__author__ = 'Joey Hwong <joeyhwong@hotmail.com>'
__url__ = 'https://github.com/ZhiweiHwang/Rhea_Chip'
__mtime__ = '20:53 01/15 2016'

from Filter import SOAPnuke
from Alignment import BWA
from GATK3 import GATK3
from Picard import Picard
from Samtools import Samtools
from MakeShell import BASH
from bamQC import DepthQC
from snpEff import snpEff
from CreateFolder import *
from CreatIndex import *
from bedAnno import *
from ExonGraph import *
from DescribeArray import DescribeArray
from qsub import Qsub
from localRun import Local
