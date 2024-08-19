#! /bin/bash

clang-format --verbose -i `find src/ -wholename "*.cpp"` `find include/ -name "*.hpp"`
