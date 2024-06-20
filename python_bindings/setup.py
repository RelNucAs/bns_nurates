#!/usr/bin/env python3

import numpy
from distutils.core import Extension, setup
from glob import glob

src_path = "../src/"

ext_list = glob(src_path + '/**/*.c', recursive=True)
ext_list.append('bns_nurates.i')
print(ext_list)
setup(
    name = "bns_nurates",
    version = "0.1",
    include_dirs = ["../include/", numpy.get_include()],
    #description = "Spin Weighted Spherical Harmonics Python Library",
    #author = "Sebastiano Bernuzzi, David Radice",
    #author_email = "sebastiano.bernuzzi@unipr.it, dradice@caltech.edu",
    ext_modules = [Extension("_bns_nurates",
      ext_list)],
#    libraries=["gsl", "gslcblas"])],
    requires = ["numpy"],
    py_modules=['bns_nurates']
    #url = "https://bitbucket.org/dradice/pyspinsph"
)
