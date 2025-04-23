#!/usr/bin/env python3

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os

BS_REAL = "double" # float

bnsnurates_mod = Extension(
    "_bnsnurates",                                                               # Name of the extension module
    swig_opts=['-python', '-c++', f'-DREAL_TYPE={BS_REAL:s}', '-I../include/'],  # SWIG options: generate Python bindings and use C++
    sources=['bnsnurates.i'],                                                    # SWIG wrapper file and C++ code
    include_dirs = ["../include/"],                                              # Include directory
    extra_compile_args=['-std=c++17', f'-DREAL_TYPE={BS_REAL:s}'],               # Additional compilation arguments
)

setup(
    name = "bnsnurates",
    version = "1.1",
    description = "Python bindings for the BNS_NURATES C++ library",
    author = "Albino Perego, David Radice, Federico Maria Guercilena, Leonardo Chiesa and Maitraya Bhattacharyya",
    author_email = "leonardo.chiesa@unitn.it",
    ext_modules = [bnsnurates_mod],
#    libraries=["gsl", "gslcblas"])],
#    requires = ["numpy"],
    py_modules=['bnsnurates'],
    url = "https://github.com/RelNucAs/bns_nurates"
)
