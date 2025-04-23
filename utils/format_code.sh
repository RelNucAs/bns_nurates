#! /bin/bash

clang-format --verbose -i ../mww.cpp `find ../include/ -name "*.hpp"` `find ../tests/ -name "*.hpp"` `find ../tests/ -name "*.cpp"`
