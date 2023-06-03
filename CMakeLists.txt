cmake_minimum_required(VERSION 3.25)
project(bns_nurates C)

set(CMAKE_C_STANDARD 11)

add_executable(bns_nurates
        main.c
        src/constants.h
        src/bns_nurates.h
        src/functions/functions.h
        src/functions/gamma.c
        src/functions/safe_exp.c
        src/integration/integration.h
        src/opacities/kernels/kernels.h
        src/opacities/weak_magnetism/weak_magnetism.h
        src/opacities/opacities.h
        src/output/outputs.h
        src/eos/eos.h
        tests/tests.h)

target_link_libraries(bns_nurates m)
