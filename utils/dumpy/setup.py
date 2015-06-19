# -*- coding: utf-8 -*-
import numpy
from numpy.distutils.core import setup, Extension

PATH_INCLUDES = [numpy.get_include()]
PATH_LIBRARIES = ['build']
LINK_LIBRARIES = []

setup(
    name = "dumpy",
    version = "1.0",
    packages = ["dumpy", "dumpy.data", "dumpy.utils", "dumpy.plots" \
                       , "dumpy.analysis"],
    )
