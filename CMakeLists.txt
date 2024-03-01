cmake_minimum_required(VERSION 3.2)
project(bns_nurates C)

set(CMAKE_C_STANDARD 11)
find_package(GSL)

add_definitions(-DSOURCE_DIR="${CMAKE_SOURCE_DIR}")

add_executable(bns_nurates
        main.c
        include/constants.h
        include/bns_nurates.h
        include/distribution.h
        src/distribution/distribution_opt_thick.c
        src/distribution/distribution_opt_thin.c
        src/distribution/distribution_total.c
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_integrals.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/lambert.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        include/opacities.h
        src/opacities/nu_abs_em_beta.c
        src/opacities/nu_scatt_iso.c
        src/output/outputs.h
        include/eos.h
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_brem.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/kernels/kernel_nes.c
        src/opacities/nu_abs_em_pair.c
        src/opacities/nu_abs_em_brem.c
        src/integration/integration_quadrature_integrate_kernels.c
        src/opacities/m1_opacities.c)

target_compile_options(bns_nurates PRIVATE -Wall -Wextra -pedantic)
target_include_directories(bns_nurates PRIVATE include)
target_include_directories(bns_nurates PRIVATE src)

target_link_libraries(bns_nurates GSL::gsl GSL::gslcblas)
target_link_libraries(bns_nurates m)

# create shared library for THC
add_library(nurates SHARED
        src/distribution/distribution_opt_thick.c
        src/distribution/distribution_opt_thin.c
        src/distribution/distribution_total.c
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/fermi_integrals.c
        src/functions/gamma.c
        src/functions/lambert.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        src/opacities/nu_abs_em_beta.c
        src/opacities/nu_scatt_iso.c
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_brem.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/kernels/kernel_nes.c
        src/opacities/nu_abs_em_pair.c
        src/opacities/nu_abs_em_brem.c
        src/integration/integration_quadrature_integrate_kernels.c
        src/opacities/m1_opacities.c)

#set_target_properties(nurates PROPERTIES PUBLIC_HEADER include/*.h)
include_directories(include)
include_directories(src)

file(GLOB HEADERS_ALL "include/*.h")
set_target_properties(nurates PROPERTIES PUBLIC_HEADER "${HEADERS_ALL}")

target_include_directories(nurates PRIVATE include)
target_include_directories(nurates PRIVATE include)
target_include_directories(nurates PRIVATE src)
target_link_libraries(nurates GSL::gsl GSL::gslcblas)
target_link_libraries(nurates m)

include(GNUInstallDirs)
install(TARGETS nurates LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nurates)

# test pair kernel T
add_executable(test_pair_t
        tests/tests_pair/test_pair_t.c
        include/constants.h
        src/functions/safe_exp.c
        include/kernels.h
        src/opacities/kernels/kernel_pair.c)

target_link_libraries(test_pair_t GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_t m)

# test pair kernel F
add_executable(test_pair_f
        tests/tests_pair/test_pair_f.c
        include/constants.h
        src/functions/safe_exp.c
        include/kernels.h
        src/opacities/kernels/kernel_pair.c)

target_link_libraries(test_pair_f GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_f m)

# test pair kernel F
add_executable(test_pair_g
        tests/tests_pair/test_pair_g.c
        include/constants.h
        src/functions/safe_exp.c
        include/kernels.h
        src/opacities/kernels/kernel_pair.c)

target_link_libraries(test_pair_g GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_g m)

# test pair kernel Psi
add_executable(test_pair_psi
        tests/tests_pair/test_pair_psi.c
        include/constants.h
        src/functions/safe_exp.c
        include/kernels.h
        src/opacities/kernels/kernel_pair.c)

target_link_libraries(test_pair_psi GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_psi m)

# test pair kernel Phi
add_executable(test_pair_phi
        tests/tests_pair/test_pair_phi.c
        include/constants.h
        src/functions/safe_exp.c
        include/kernels.h
        src/opacities/kernels/kernel_pair.c
)

target_link_libraries(test_pair_phi GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_phi m)

# test pair kernel Kernels
add_executable(test_pair_kernels
        tests/tests_pair/test_pair_kernels.c
        include/constants.h
        src/functions/safe_exp.c
        include/kernels.h
        src/opacities/kernels/kernel_pair.c
)

target_link_libraries(test_pair_kernels GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_kernels m)

# test pair kernel Kernels phi integrated
add_executable(test_pair_kernels_phi_mu_integrated
        tests/tests_pair/test_pair_kernels_phi_mu_integrated.c
        include/constants.h
        src/functions/safe_exp.c
        include/kernels.h
        src/opacities/kernels/kernel_pair.c
)

target_link_libraries(test_pair_kernels_phi_mu_integrated GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_kernels_phi_mu_integrated m)

# test pair kernel opacities
add_executable(test_pair_opacities
        tests/tests_pair/test_pair_opacities.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        include/opacities.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/fermi_integrals.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        src/opacities/kernels/kernel_pair.c
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/nu_abs_em_pair.c
        src/integration/integration_quadrature_integrate_kernels.c
)

target_link_libraries(test_pair_opacities GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_opacities m)

# test beta processes spectral rates
add_executable(test_beta
        tests/tests_beta/test_beta.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/opacities.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/nu_abs_em_pair.c
        src/integration/integration_quadrature_integrate_kernels.c
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        src/opacities/nu_abs_em_beta.c
)

target_link_libraries(test_beta GSL::gsl GSL::gslcblas)
target_link_libraries(test_beta m)

# test isoenergetic scattering kernel
add_executable(test_iso
        tests/tests_iso/test_iso.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/opacities.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/nu_abs_em_pair.c
        src/integration/integration_quadrature_integrate_kernels.c
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        src/opacities/nu_scatt_iso.c
)

target_link_libraries(test_iso GSL::gsl GSL::gslcblas)
target_link_libraries(test_iso m)

# test NN bremsstrahlung kernel
add_executable(test_brem
        tests/tests_brem/test_brem.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/opacities.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/nu_abs_em_pair.c
        src/integration/integration_quadrature_integrate_kernels.c
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        src/opacities/nu_abs_em_brem.c
        src/opacities/kernels/kernel_brem.c
)

target_link_libraries(test_brem GSL::gsl GSL::gslcblas)
target_link_libraries(test_brem m)

# test gauss legendre quadrature 1d
add_executable(test_gauss_legendre_1d
        tests/tests_quadrature_integrate/test_gauss_legendre_1d.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        include/integration.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
)

target_link_libraries(test_gauss_legendre_1d GSL::gsl GSL::gslcblas)
target_link_libraries(test_gauss_legendre_1d m)

# test gauss legendre quadrature 2d
add_executable(test_gauss_legendre_2d
        tests/tests_quadrature_integrate/test_gauss_legendre_2d.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        include/integration.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
)

target_link_libraries(test_gauss_legendre_2d GSL::gsl GSL::gslcblas)
target_link_libraries(test_gauss_legendre_2d m)

# test gauss legendre quadrature 3d
add_executable(test_gauss_legendre_3d
        tests/tests_quadrature_integrate/test_gauss_legendre_3d.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        include/integration.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
)

target_link_libraries(test_gauss_legendre_3d GSL::gsl GSL::gslcblas)
target_link_libraries(test_gauss_legendre_3d m)

# test gauss laguerre quadrature 1d
add_executable(test_gauss_laguerre_1d
        tests/tests_quadrature_integrate/test_gauss_laguerre_1d.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        include/integration.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
)

target_link_libraries(test_gauss_laguerre_1d GSL::gsl GSL::gslcblas)
target_link_libraries(test_gauss_laguerre_1d m)

# test quadrature I/O
add_executable(test_quadrature_io
        tests/tests_quadrature_integrate/test_quadrature_io.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        include/integration.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
)

target_link_libraries(test_quadrature_io GSL::gsl GSL::gslcblas)
target_link_libraries(test_quadrature_io m)

# test 1d integration
add_executable(test_integration_1d
        tests/tests_quadrature_integrate/test_integration_1d.c
        tests/tests_quadrature_integrate/gsl_routines/integration_gsl.h
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        include/integration.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
)

target_link_libraries(test_integration_1d GSL::gsl GSL::gslcblas)
target_link_libraries(test_integration_1d m)

# test 2d integration
add_executable(test_integration_2d
        tests/tests_quadrature_integrate/test_integration_2d.c
        tests/tests_quadrature_integrate/gsl_routines/integration_gsl.h
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        include/integration.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
)

target_link_libraries(test_integration_2d GSL::gsl GSL::gslcblas)
target_link_libraries(test_integration_2d m)

# test 2d integration
add_executable(test_integration_3d
        tests/tests_quadrature_integrate/test_integration_3d.c
        tests/tests_quadrature_integrate/gsl_routines/integration_gsl.h
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        include/integration.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
)

target_link_libraries(test_integration_3d GSL::gsl GSL::gslcblas)
target_link_libraries(test_integration_3d m)

# test integration special
add_executable(test_integration_special
        tests/tests_quadrature_integrate/test_integration_special.c
        tests/tests_quadrature_integrate/gsl_routines/integration_gsl.h
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/integration/integration_quadrature_integrate_kernels.c
)

target_link_libraries(test_integration_special GSL::gsl GSL::gslcblas)
target_link_libraries(test_integration_special m)

# test integration special 1d
add_executable(test_integration_multi_func_1d
        tests/tests_quadrature_integrate/test_integration_multi_func_1d.c
        tests/tests_quadrature_integrate/gsl_routines/integration_gsl.h
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/integration/integration_quadrature_integrate_kernels.c
)

target_link_libraries(test_integration_multi_func_1d GSL::gsl GSL::gslcblas)
target_link_libraries(test_integration_multi_func_1d m)

# test integration special 2d
add_executable(test_integration_multi_func_2d
        tests/tests_quadrature_integrate/test_integration_multi_func_2d.c
        tests/tests_quadrature_integrate/gsl_routines/integration_gsl.h
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        include/constants.h
        include/bns_nurates.h
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/integration/integration_quadrature_integrate_kernels.c
)

target_link_libraries(test_integration_multi_func_2d GSL::gsl GSL::gslcblas)
target_link_libraries(test_integration_multi_func_2d m)

# test opacities electron capture
add_executable(test_opacities_electron_capture
        include/constants.h
        include/bns_nurates.h
        include/distribution.h
        src/distribution/distribution_opt_thick.c
        src/distribution/distribution_opt_thin.c
        src/distribution/distribution_total.c
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_integrals.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/lambert.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        include/opacities.h
        src/opacities/nu_abs_em_beta.c
        src/opacities/nu_scatt_iso.c
        src/output/outputs.h
        include/eos.h
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_brem.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/kernels/kernel_nes.c
        src/opacities/nu_abs_em_pair.c
        src/opacities/nu_abs_em_brem.c
        src/integration/integration_quadrature_integrate_kernels.c
        src/opacities/m1_opacities.c
        tests/tests_opacities_m1/test_opacities_electron_capture.c
        tests/tests_opacities_m1/test_m1_opacities.c
        tests/tests.h 
)

target_link_libraries(test_opacities_electron_capture GSL::gsl GSL::gslcblas)
target_link_libraries(test_opacities_electron_capture m)

# test opacities pair
add_executable(test_opacities_pair
        include/constants.h
        include/bns_nurates.h
        include/distribution.h
        src/distribution/distribution_opt_thick.c
        src/distribution/distribution_opt_thin.c
        src/distribution/distribution_total.c
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_integrals.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/lambert.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        include/opacities.h
        src/opacities/nu_abs_em_beta.c
        src/opacities/nu_scatt_iso.c
        src/output/outputs.h
        include/eos.h
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_brem.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/kernels/kernel_nes.c
        src/opacities/nu_abs_em_pair.c
        src/opacities/nu_abs_em_brem.c
        src/integration/integration_quadrature_integrate_kernels.c
        src/opacities/m1_opacities.c
        tests/tests_opacities_m1/test_opacities_pair.c
        tests/tests_opacities_m1/test_m1_opacities.c
        tests/tests.h 
)

target_link_libraries(test_opacities_pair GSL::gsl GSL::gslcblas)
target_link_libraries(test_opacities_pair m)

# test opacities bremsstrahlung
add_executable(test_opacities_brem
        include/constants.h
        include/bns_nurates.h
        include/distribution.h
        src/distribution/distribution_opt_thick.c
        src/distribution/distribution_opt_thin.c
        src/distribution/distribution_total.c
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_integrals.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/lambert.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        include/opacities.h
        src/opacities/nu_abs_em_beta.c
        src/opacities/nu_scatt_iso.c
        src/output/outputs.h
        include/eos.h
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_brem.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/kernels/kernel_nes.c
        src/opacities/nu_abs_em_pair.c
        src/opacities/nu_abs_em_brem.c
        src/integration/integration_quadrature_integrate_kernels.c
        src/opacities/m1_opacities.c
        tests/tests_opacities_m1/test_opacities_brem.c
        tests/tests_opacities_m1/test_m1_opacities.c
        tests/tests.h
)

target_link_libraries(test_opacities_brem GSL::gsl GSL::gslcblas)
target_link_libraries(test_opacities_brem m)

# test opacities inelastic scattering
add_executable(test_opacities_inelastic_scatt
        include/constants.h
        include/bns_nurates.h
        include/distribution.h
        src/distribution/distribution_opt_thick.c
        src/distribution/distribution_opt_thin.c
        src/distribution/distribution_total.c
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_integrals.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/lambert.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        include/opacities.h
        src/opacities/nu_abs_em_beta.c
        src/opacities/nu_scatt_iso.c
        src/output/outputs.h
        include/eos.h
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_brem.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/kernels/kernel_nes.c
        src/opacities/nu_abs_em_pair.c
        src/opacities/nu_abs_em_brem.c
        src/integration/integration_quadrature_integrate_kernels.c
        src/opacities/m1_opacities.c
        tests/tests_opacities_m1/test_opacities_inelastic_scatt.c
        tests/tests_opacities_m1/test_m1_opacities.c
        tests/tests.h 
)

target_link_libraries(test_opacities_inelastic_scatt GSL::gsl GSL::gslcblas)
target_link_libraries(test_opacities_inelastic_scatt m)

# test opacities isoscatt
add_executable(test_opacities_isoscatt
        include/constants.h
        include/bns_nurates.h
        include/distribution.h
        src/distribution/distribution_opt_thick.c
        src/distribution/distribution_opt_thin.c
        src/distribution/distribution_total.c
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_integrals.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/lambert.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        include/opacities.h
        src/opacities/nu_abs_em_beta.c
        src/opacities/nu_scatt_iso.c
        src/output/outputs.h
        include/eos.h
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_brem.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/kernels/kernel_nes.c
        src/opacities/nu_abs_em_pair.c
        src/opacities/nu_abs_em_brem.c
        src/integration/integration_quadrature_integrate_kernels.c
        src/opacities/m1_opacities.c
        tests/tests_opacities_m1/test_opacities_isoscatt.c
        tests/tests_opacities_m1/test_m1_opacities.c
        tests/tests.h
)

target_link_libraries(test_opacities_isoscatt GSL::gsl GSL::gslcblas)
target_link_libraries(test_opacities_isoscatt m)

# test opacities all reactions
add_executable(test_opacities_all_reactions
        include/constants.h
        include/bns_nurates.h
        include/distribution.h
        src/distribution/distribution_opt_thick.c
        src/distribution/distribution_opt_thin.c
        src/distribution/distribution_total.c
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_integrals.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/lambert.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        include/opacities.h
        src/opacities/nu_abs_em_beta.c
        src/opacities/nu_scatt_iso.c
        src/output/outputs.h
        include/eos.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_brem.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/kernels/kernel_nes.c
        src/opacities/nu_abs_em_pair.c
        src/opacities/nu_abs_em_brem.c
        src/integration/integration_quadrature_integrate_kernels.c
        src/opacities/m1_opacities.c
        tests/tests_opacities_m1/test_opacities_all_reactions.c
        tests/tests_opacities_m1/test_m1_opacities.c
        tests/tests.h 
)
target_link_libraries(test_opacities_all_reactions GSL::gsl GSL::gslcblas)
target_link_libraries(test_opacities_all_reactions m)

target_link_libraries(test_opacities_all_reactions GSL::gsl GSL::gslcblas)
target_link_libraries(test_opacities_all_reactions m)

add_executable(test_opacities_thc
        include/constants.h
        include/bns_nurates.h
        include/distribution.h
        src/distribution/distribution_opt_thick.c
        src/distribution/distribution_opt_thin.c
        src/distribution/distribution_total.c
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_integrals.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/lambert.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        include/opacities.h
        src/opacities/nu_abs_em_beta.c
        src/opacities/nu_scatt_iso.c
        src/output/outputs.h
        include/eos.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_brem.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/kernels/kernel_nes.c
        src/opacities/nu_abs_em_pair.c
        src/opacities/nu_abs_em_brem.c
        src/integration/integration_quadrature_integrate_kernels.c
        src/opacities/m1_opacities.c
        tests/tests_opacities_m1/test_opacities_thc.c
        tests/tests_opacities_m1/test_m1_opacities.c
        tests/tests.h
)
target_link_libraries(test_opacities_thc GSL::gsl GSL::gslcblas)
target_link_libraries(test_opacities_thc m)

target_link_libraries(test_opacities_thc GSL::gsl GSL::gslcblas)
target_link_libraries(test_opacities_thc m)

# test m1 distribution
add_executable(test_m1_distribution
        include/constants.h
        include/bns_nurates.h
        include/distribution.h
        src/distribution/distribution_opt_thick.c
        src/distribution/distribution_opt_thin.c
        src/distribution/distribution_total.c
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_integrals.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/lambert.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        include/opacities.h
        src/opacities/nu_abs_em_beta.c
        src/opacities/nu_scatt_iso.c
        src/output/outputs.h
        include/eos.h
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.c
        tests/tests_quadrature_integrate/gsl_routines/integration_GSL.c
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_brem.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/kernels/kernel_nes.c
        src/opacities/nu_abs_em_pair.c
        src/opacities/nu_abs_em_brem.c
        src/integration/integration_quadrature_integrate_kernels.c
        src/opacities/m1_opacities.c
        tests/tests_distribution/test_m1_distribution.c
)

target_link_libraries(test_m1_distribution GSL::gsl GSL::gslcblas)
target_link_libraries(test_m1_distribution m)

# test quadrature convergence spectral
add_executable(test_quadrature_convergence_spectral
        include/constants.h
        include/bns_nurates.h
        include/distribution.h
        src/distribution/distribution_opt_thick.c
        src/distribution/distribution_opt_thin.c
        src/distribution/distribution_total.c
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_integrals.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/lambert.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        include/opacities.h
        src/opacities/nu_abs_em_beta.c
        src/opacities/nu_scatt_iso.c
        src/output/outputs.h
        include/eos.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_brem.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/kernels/kernel_nes.c
        src/opacities/nu_abs_em_pair.c
        src/opacities/nu_abs_em_brem.c
        src/integration/integration_quadrature_integrate_kernels.c
        src/opacities/m1_opacities.c
        tests/tests_quadrature_convergence/test_quadrature_convergence_spectral.c
        tests/tests.h 
)

target_link_libraries(test_quadrature_convergence_spectral m)
target_link_libraries(test_quadrature_convergence_spectral GSL::gsl GSL::gslcblas)

# test quadrature convergence integrated
add_executable(test_quadrature_convergence_integrated
        include/constants.h
        include/bns_nurates.h
        include/distribution.h
        src/distribution/distribution_opt_thick.c
        src/distribution/distribution_opt_thin.c
        src/distribution/distribution_total.c
        include/functions.h
        src/functions/digamma.c
        src/functions/fermi_integrals.c
        src/functions/fermi_distr.c
        src/functions/gamma.c
        src/functions/lambert.c
        src/functions/mnewt.c
        src/functions/safe_exp.c
        src/functions/theta.c
        src/functions/m1matrix.c
        include/integration.h
        include/kernels.h
        include/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
        include/opacities.h
        src/opacities/nu_abs_em_beta.c
        src/opacities/nu_scatt_iso.c
        src/output/outputs.h
        include/eos.h
        src/integration/integration_quadrature_schemes.c
        src/integration/integration_quadrature_io.c
        src/integration/integration_quadrature_integrate.c
        src/opacities/kernels/kernel_brem.c
        src/opacities/kernels/kernel_pair.c
        src/opacities/kernels/kernel_nes.c
        src/opacities/nu_abs_em_pair.c
        src/opacities/nu_abs_em_brem.c
        src/integration/integration_quadrature_integrate_kernels.c
        src/opacities/m1_opacities.c
        tests/tests_quadrature_convergence/test_quadrature_convergence_integrated.c
        tests/tests.h 
)

target_link_libraries(test_quadrature_convergence_integrated m)
target_link_libraries(test_quadrature_convergence_integrated GSL::gsl GSL::gslcblas)

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
add_test(NAME test_pair_kernels_phi_mu_integrated
        COMMAND $<TARGET_FILE:test_pair_kernels_phi_mu_integrated>)
add_test(NAME test_pair_opacities
        COMMAND $<TARGET_FILE:test_pair_opacities>)
add_test(NAME test_gauss_legendre_1d
        COMMAND $<TARGET_FILE:test_gauss_legendre_1d>)
add_test(NAME test_gauss_legendre_2d
        COMMAND $<TARGET_FILE:test_gauss_legendre_2d>)
add_test(NAME test_gauss_legendre_3d
        COMMAND $<TARGET_FILE:test_gauss_legendre_3d>)
add_test(NAME test_gauss_laguerre_1d
        COMMAND $<TARGET_FILE:test_gauss_laguerre_1d>)
add_test(NAME test_quadrature_io
        COMMAND $<TARGET_FILE:test_quadrature_io>)
add_test(NAME test_integration_1d
        COMMAND $<TARGET_FILE:test_integration_1d>)
add_test(NAME test_integration_2d
        COMMAND $<TARGET_FILE:test_integration_2d>)
add_test(NAME test_integration_3d
        COMMAND $<TARGET_FILE:test_integration_3d>)
add_test(NAME test_integration_special
        COMMAND $<TARGET_FILE:test_integration_special>)
add_test(NAME test_integration_multi_func_1d
        COMMAND $<TARGET_FILE:test_integration_multi_func_1d>)
add_test(NAME test_integration_multi_func_2d
        COMMAND $<TARGET_FILE:test_integration_multi_func_2d>)
add_test(NAME test_opacities_electron_capture
        COMMAND $<TARGET_FILE:test_opacities_electron_capture>)
add_test(NAME test_opacities_pair
        COMMAND $<TARGET_FILE:test_opacities_pair>)
add_test(NAME test_opacities_brem
        COMMAND $<TARGET_FILE:test_opacities_brem>)
add_test(NAME test_opacities_isoscatt
        COMMAND $<TARGET_FILE:test_opacities_isoscatt>)
add_test(NAME test_m1_distribution
        COMMAND $<TARGET_FILE:test_m1_distribution>)
add_test(NAME test_opacities_thc
        COMMAND $<TARGET_FILE:test_opacities_thc>)