cmake_minimum_required(VERSION 3.25)
project(bns_nurates C)

set(CMAKE_C_STANDARD 11)
find_package(GSL)

add_executable(bns_nurates
        main.c
        src/constants.h
        src/bns_nurates.h
        src/functions/functions.h
        src/functions/digamma.c
        src/functions/fermi.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/integration/integration.h
        src/opacities/kernels/kernels.h
        src/opacities/kernels/kernel_iso_scatt.c
        src/opacities/weak_magnetism/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        src/opacities/opacities.h
        src/opacities/nu_abs_em_nuclabs.c
        src/output/outputs.h
        src/eos/eos.h
        tests/tests.h
        tests/GSL_code/integration_GSL.c
        tests/GSL_code/mnewt_GSL.c
        tests/GSL_code/integration_GSL.c
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        tests/tests_quadrature.c
        src/opacities/kernels/kernel_brem.c
        src/opacities/kernels/kernel_pair.c
        tests/tests_integration.c
        src/opacities/nu_abs_em_pair.c
        src/opacities/nu_abs_em_brem.c
        src/integration/integration_quadrature_integrate_kernels.c tests/test_pair.c tests/test_brem.c)

target_link_libraries(bns_nurates GSL::gsl GSL::gslcblas)
target_link_libraries(bns_nurates m)

# test pair kernel T
add_executable(test_pair_t
        tests/tests_pair/test_pair_t.c
        src/constants.h
        src/opacities/kernels/kernels.h
        src/opacities/kernels/kernel_pair.c)

target_link_libraries(test_pair_t GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_t m)

# test pair kernel F
add_executable(test_pair_f
        tests/tests_pair/test_pair_f.c
        src/constants.h
        src/opacities/kernels/kernels.h
        src/opacities/kernels/kernel_pair.c)

target_link_libraries(test_pair_f GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_f m)

# test pair kernel F
add_executable(test_pair_g
        tests/tests_pair/test_pair_g.c
        src/constants.h
        src/opacities/kernels/kernels.h
        src/opacities/kernels/kernel_pair.c)

target_link_libraries(test_pair_g GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_g m)

# test pair kernel Psi
add_executable(test_pair_psi
        tests/tests_pair/test_pair_psi.c
        src/constants.h
        src/opacities/kernels/kernels.h
        src/opacities/kernels/kernel_pair.c)

target_link_libraries(test_pair_psi GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_psi m)

# test pair kernel Phi
add_executable(test_pair_phi
        tests/tests_pair/test_pair_phi.c
        src/constants.h
        src/opacities/kernels/kernels.h
        src/opacities/kernels/kernel_pair.c
        )

target_link_libraries(test_pair_phi GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_phi m)

# test pair kernel Kernels
add_executable(test_pair_kernels
        tests/tests_pair/test_pair_kernels.c
        src/constants.h
        src/opacities/kernels/kernels.h
        src/opacities/kernels/kernel_pair.c
        )

target_link_libraries(test_pair_kernels GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_kernels m)

# test pair kernel opacities
add_executable(test_pair_opacities
        tests/tests_pair/test_pair_opacities.c
        src/constants.h
        src/bns_nurates.h
        src/functions/functions.h
        src/functions/digamma.c
        src/functions/fermi.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/integration/integration.h
        src/opacities/kernels/kernels.h
        src/opacities/opacities.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/nu_abs_em_pair.c
        src/integration/integration_quadrature_integrate_kernels.c
        )

target_link_libraries(test_pair_opacities GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_opacities m)

enable_testing()

# unit tests
add_test(NAME test_pair_t
        COMMAND $<TARGET_FILE:test_pair_t>)
add_test(NAME test_pair_f
        COMMAND $<TARGET_FILE:test_pair_f>)
add_test(NAME test_pair_g
        COMMAND $<TARGET_FILE:test_pair_g>)
add_test(NAME test_pair_psi
        COMMAND $<TARGET_FILE:test_pair_psi>)
add_test(NAME test_pair_phi
        COMMAND $<TARGET_FILE:test_pair_phi>)
add_test(NAME test_pair_kernels
        COMMAND $<TARGET_FILE:test_pair_kernels>)
add_test(NAME test_pair_opacities
        COMMAND $<TARGET_FILE:test_pair_opacities>)