#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

import os
from collections import defaultdict

class JavaTool(object):
	def __init__(self, jar_path, tmp_path=os.getcwd(), **kwargs):
		kwargs = defaultdict(lambda :None, kwargs)
		self.tools = os.path.abspath(kwargs["java_path"]) if kwargs["java_path"] else "java"
		self.options = self.tools
		self.options += " -Xmx%s" % kwargs["max_memory"] if kwargs["max_memory"] else ""
		self.options += " -Djava.io.tmpdir=%s -jar %s" % (os.path.abspath(tmp_path), os.path.abspath(jar_path))
