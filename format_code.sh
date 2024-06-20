#! /bin/bash

clang-format --verbose -i `find src/ -wholename "*.c"` `find include/ -name "*.h"`
