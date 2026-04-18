#! /bin/bash

clang-format --verbose -i ../mwe.cpp `find ../include/ -name "*.hpp"` `find ../tests/ -name "*.hpp"` `find ../tests/ -name "*.cpp"`
