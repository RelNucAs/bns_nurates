#! /bin/bash

clang-format --verbose -i `find ../include/ -name "*.hpp"` `find ../tests/ -name "*.hpp"` `find ../tests/ -name "*.cpp"`
