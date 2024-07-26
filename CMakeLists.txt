cmake_minimum_required(VERSION 3.2)
project(bns_nurates VERSION 0.1 LANGUAGES C CXX)

# Kokkos require C++17 standard as minimum
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Force flag to denote C++ standard, as deprecated c++1z is not supported by nvcc.
set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++17")

# Kokkos recommends NOT using default compiler extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# following policy required to enable target_link_libraries() below
cmake_policy(SET CMP0079 NEW)

find_package(GSL)

# set OpenMP macro (true/false)
set(ENABLE_OPENMP OFF)
find_package(OpenMP COMPONENTS CXX)
if (NOT OpenMP_FOUND)
  message(FATAL_ERROR "OpenMP package is required but could not be found.")
 endif()
set(ENABLE_OPENMP ON)
set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "Compile Kokkos with OpenMP")
if (ENABLE_OPENMP)
  set(OPENMP_PARALLEL_ENABLED 1)
else()
  set(OPENMP_PARALLEL_ENABLED 0)
endif()

if (${GSL_FOUND})
	add_compile_definitions(GSL_INCLUDES_H_)
endif()

add_definitions(-DSOURCE_DIR="${CMAKE_SOURCE_DIR}")

add_executable(bns_nurates
        main.c
		include/constants.hpp
		include/bns_nurates.hpp
		include/distribution.hpp
		include/functions.hpp
		include/integration.hpp
		include/kernels.hpp
		include/weak_magnetism.hpp
		include/opacities.hpp
		src/output/outputs.hpp
		include/eos.hpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/kernel_nes.hpp
		include/kernel_brem.hpp)

add_subdirectory(kokkos)
include_directories(${Kokkos_INCLUDE_DIRS_RET})
target_link_libraries(bns_nurates Kokkos::kokkos)
target_link_libraries(bns_nurates OpenMP::OpenMP_CXX)

#target_compile_options(bns_nurates PRIVATE -Wall -Wextra)
target_include_directories(bns_nurates PRIVATE include)
target_include_directories(bns_nurates PRIVATE src)

if (${GSL_FOUND})
	target_link_libraries(bns_nurates GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(bns_nurates m)

# create shared library for THC
add_library(nurates SHARED
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/kernel_nes.hpp)

#set_target_properties(nurates PROPERTIES PUBLIC_HEADER include/*.h)
include_directories(include)
include_directories(src)

file(GLOB HEADERS_ALL "include/*.h")
set_target_properties(nurates PROPERTIES PUBLIC_HEADER "${HEADERS_ALL}")

target_include_directories(nurates PRIVATE include)
target_include_directories(nurates PRIVATE include)
target_include_directories(nurates PRIVATE src)
if (${GSL_FOUND})
	target_link_libraries(nurates GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(nurates m)
target_link_libraries(nurates Kokkos::kokkos)
target_link_libraries(nurates OpenMP::OpenMP_CXX)

include(GNUInstallDirs)
install(TARGETS nurates LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nurates)

# test pair kernel T
add_executable(test_pair_t
		tests/tests_pair/test_pair_t.cpp
		include/constants.hpp
		include/kernels.hpp
		include/kernel_pair.hpp)
target_link_libraries(test_pair_t Kokkos::kokkos)
target_link_libraries(test_pair_t OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_pair_t GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_pair_t m)

# test pair kernel F
add_executable(test_pair_f
		tests/tests_pair/test_pair_f.cpp
		include/constants.hpp
		include/kernels.hpp
		include/kernel_pair.hpp)

target_link_libraries(test_pair_f Kokkos::kokkos)
target_link_libraries(test_pair_f OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_pair_f GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_pair_f m)

# test pair kernel F
add_executable(test_pair_g
		tests/tests_pair/test_pair_g.cpp
		include/constants.hpp
		include/kernels.hpp
		include/kernel_pair.hpp)
target_link_libraries(test_pair_g Kokkos::kokkos)
target_link_libraries(test_pair_g OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_pair_g GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_pair_g m)

# test pair kernel Psi
add_executable(test_pair_psi
		tests/tests_pair/test_pair_psi.cpp
		include/constants.hpp
		include/kernels.hpp
		include/kernel_pair.hpp)
target_link_libraries(test_pair_psi Kokkos::kokkos)
target_link_libraries(test_pair_psi OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_pair_psi GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_pair_psi m)

# test pair kernel Phi
add_executable(test_pair_phi
		tests/tests_pair/test_pair_phi.cpp
		include/constants.hpp
		include/kernels.hpp
		include/kernel_pair.hpp
)
target_link_libraries(test_pair_phi Kokkos::kokkos)
target_link_libraries(test_pair_phi OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_pair_phi GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_pair_phi m)

if (${GSL_FOUND})
# test pair kernel Kernels
add_executable(test_pair_kernels
		tests/tests_pair/test_pair_kernels.cpp
		include/constants.hpp
		include/kernels.hpp
		include/kernel_pair.hpp
)
target_link_libraries(test_pair_kernels Kokkos::kokkos)
target_link_libraries(test_pair_kernels OpenMP::OpenMP_CXX)
target_link_libraries(test_pair_kernels GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_kernels m)
endif()


if (${GSL_FOUND})
# test pair kernel Kernels phi integrated
add_executable(test_pair_kernels_phi_mu_integrated
		tests/tests_pair/test_pair_kernels_phi_mu_integrated.cpp
		include/constants.hpp
		include/kernels.hpp
		include/kernel_pair.hpp
)
target_link_libraries(test_pair_kernels_phi_mu_integrated Kokkos::kokkos)
target_link_libraries(test_pair_kernels_phi_mu_integrated OpenMP::OpenMP_CXX)
target_link_libraries(test_pair_kernels_phi_mu_integrated GSL::gsl GSL::gslcblas)
target_link_libraries(test_pair_kernels_phi_mu_integrated m)
endif()


# test pair kernel opacities
add_executable(test_pair_opacities
		tests/tests_pair/test_pair_opacities.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/opacities.hpp
		include/integration.hpp
		include/kernels.hpp
		include/kernel_pair.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
)
target_link_libraries(test_pair_opacities Kokkos::kokkos)
target_link_libraries(test_pair_opacities OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_pair_opacities GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_pair_opacities m)

# test beta processes spectral rates
add_executable(test_beta
		tests/tests_beta/test_beta.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/kernels.hpp
		include/weak_magnetism.hpp
)

target_link_libraries(test_beta Kokkos::kokkos)
target_link_libraries(test_beta OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_beta GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_beta m)

# test isoenergetic scattering kernel
add_executable(test_iso
		tests/tests_iso/test_iso.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/integration.hpp
		include/kernels.hpp
		include/opacities.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/weak_magnetism.hpp
)

target_link_libraries(test_iso Kokkos::kokkos)
target_link_libraries(test_iso OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_iso GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_iso m)

# test NN bremsstrahlung kernel
add_executable(test_brem
		tests/tests_brem/test_brem.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/integration.hpp
		include/kernels.hpp
		include/opacities.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/weak_magnetism.hpp
)

target_link_libraries(test_brem Kokkos::kokkos)
target_link_libraries(test_brem OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_brem GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_brem m)

# test gauss legendre quadrature 1d
add_executable(test_gauss_legendre_1d
		tests/tests_quadrature_integrate/test_gauss_legendre_1d.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/integration.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
)

target_link_libraries(test_gauss_legendre_1d Kokkos::kokkos)
target_link_libraries(test_gauss_legendre_1d OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_gauss_legendre_1d GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_gauss_legendre_1d m)

# test gauss legendre quadrature 2d
add_executable(test_gauss_legendre_2d
		tests/tests_quadrature_integrate/test_gauss_legendre_2d.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/integration.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
)

target_link_libraries(test_gauss_legendre_2d Kokkos::kokkos)
target_link_libraries(test_gauss_legendre_2d OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_gauss_legendre_2d GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_gauss_legendre_2d m)

# test gauss legendre quadrature 3d
add_executable(test_gauss_legendre_3d
		tests/tests_quadrature_integrate/test_gauss_legendre_3d.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/integration.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
)

target_link_libraries(test_gauss_legendre_3d Kokkos::kokkos)
target_link_libraries(test_gauss_legendre_3d OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_gauss_legendre_3d GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_gauss_legendre_3d m)

# test gauss laguerre quadrature 1d
add_executable(test_gauss_laguerre_1d
		tests/tests_quadrature_integrate/test_gauss_laguerre_1d.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/integration.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
)

target_link_libraries(test_gauss_laguerre_1d Kokkos::kokkos)
target_link_libraries(test_gauss_laguerre_1d OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_gauss_laguerre_1d GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_gauss_laguerre_1d m)

# test quadrature I/O
add_executable(test_quadrature_io
		tests/tests_quadrature_integrate/test_quadrature_io.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/integration.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
)
target_link_libraries(test_quadrature_io Kokkos::kokkos)
target_link_libraries(test_quadrature_io OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_quadrature_io GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_quadrature_io m)

# test 1d integration
add_executable(test_integration_1d
		tests/tests_quadrature_integrate/test_integration_1d.cpp
		tests/tests_quadrature_integrate/gsl_routines/integration_gsl.hpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/integration.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
)

target_link_libraries(test_integration_1d Kokkos::kokkos)
target_link_libraries(test_integration_1d OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_integration_1d GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_integration_1d m)

# test 2d integration
add_executable(test_integration_2d
		tests/tests_quadrature_integrate/test_integration_2d.cpp
		tests/tests_quadrature_integrate/gsl_routines/integration_gsl.hpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/integration.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
)

target_link_libraries(test_integration_2d Kokkos::kokkos)
target_link_libraries(test_integration_2d OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_integration_2d GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_integration_2d m)

# test 2d integration
add_executable(test_integration_3d
		tests/tests_quadrature_integrate/test_integration_3d.cpp
		tests/tests_quadrature_integrate/gsl_routines/integration_gsl.hpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/integration.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
)

target_link_libraries(test_integration_3d Kokkos::kokkos)
target_link_libraries(test_integration_3d OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_integration_3d GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_integration_3d m)

# test integration special
add_executable(test_integration_special
		tests/tests_quadrature_integrate/test_integration_special.cpp
		tests/tests_quadrature_integrate/gsl_routines/integration_gsl.hpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/integration.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
)

target_link_libraries(test_integration_special Kokkos::kokkos)
target_link_libraries(test_integration_special OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_integration_special GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_integration_special m)

# test integration special 1d
add_executable(test_integration_multi_func_1d
		tests/tests_quadrature_integrate/test_integration_multi_func_1d.cpp
		tests/tests_quadrature_integrate/gsl_routines/integration_gsl.hpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/integration.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
)

target_link_libraries(test_integration_multi_func_1d Kokkos::kokkos)
target_link_libraries(test_integration_multi_func_1d OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_integration_multi_func_1d GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_integration_multi_func_1d m)

# test integration special 2d
add_executable(test_integration_multi_func_2d
		tests/tests_quadrature_integrate/test_integration_multi_func_2d.cpp
		tests/tests_quadrature_integrate/gsl_routines/integration_gsl.hpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		include/constants.hpp
		include/bns_nurates.hpp
		include/functions.hpp
		include/integration.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
)

target_link_libraries(test_integration_multi_func_2d Kokkos::kokkos)
target_link_libraries(test_integration_multi_func_2d OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_integration_multi_func_2d GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_integration_multi_func_2d m)

# test opacities electron capture
add_executable(test_opacities_electron_capture
		include/constants.hpp
		include/bns_nurates.hpp
		include/distribution.hpp
		include/functions.hpp
		include/integration.hpp
		include/kernels.hpp
		include/weak_magnetism.hpp
		include/opacities.hpp
		src/output/outputs.hpp
		include/eos.hpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/kernel_nes.hpp
		tests/tests_opacities_m1/test_opacities_electron_capture.cpp
        tests/tests_opacities_m1/test_m1_opacities.cpp
		tests/tests.hpp
)
target_link_libraries(test_opacities_electron_capture Kokkos::kokkos)
target_link_libraries(test_opacities_electron_capture OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_opacities_electron_capture GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_opacities_electron_capture m)

# test opacities pair
add_executable(test_opacities_pair
		include/constants.hpp
		include/bns_nurates.hpp
		include/distribution.hpp
		include/functions.hpp
		include/integration.hpp
		include/kernels.hpp
		include/weak_magnetism.hpp
		include/opacities.hpp
		src/output/outputs.hpp
		include/eos.hpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/kernel_nes.hpp
		tests/tests_opacities_m1/test_opacities_pair.cpp
        tests/tests_opacities_m1/test_m1_opacities.cpp
		tests/tests.hpp
)
target_link_libraries(test_opacities_pair Kokkos::kokkos)
target_link_libraries(test_opacities_pair OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_opacities_pair GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_opacities_pair m)

# test opacities bremsstrahlung
add_executable(test_opacities_brem
		include/constants.hpp
		include/bns_nurates.hpp
		include/distribution.hpp
		include/functions.hpp
		include/integration.hpp
		include/kernels.hpp
		include/weak_magnetism.hpp
		include/opacities.hpp
		src/output/outputs.hpp
		include/eos.hpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/kernel_nes.hpp
		tests/tests_opacities_m1/test_opacities_brem.cpp
        tests/tests_opacities_m1/test_m1_opacities.cpp
		tests/tests.hpp
		include/kernel_brem.hpp
)

target_link_libraries(test_opacities_brem Kokkos::kokkos)
target_link_libraries(test_opacities_brem OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_opacities_brem GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_opacities_brem m)

# test opacities inelastic scattering
add_executable(test_opacities_inelastic_scatt
		include/constants.hpp
		include/bns_nurates.hpp
		include/distribution.hpp
		include/functions.hpp
		include/integration.hpp
		include/kernels.hpp
		include/weak_magnetism.hpp
		include/opacities.hpp
		src/output/outputs.hpp
		include/eos.hpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/kernel_nes.hpp
		tests/tests_opacities_m1/test_opacities_inelastic_scatt.cpp
        tests/tests_opacities_m1/test_m1_opacities.cpp
		tests/tests.hpp
)
target_link_libraries(test_opacities_inelastic_scatt Kokkos::kokkos)
target_link_libraries(test_opacities_inelastic_scatt OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_opacities_inelastic_scatt GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_opacities_inelastic_scatt m)

# test opacities isoscatt
add_executable(test_opacities_isoscatt
		include/constants.hpp
		include/bns_nurates.hpp
		include/distribution.hpp
		include/functions.hpp
		include/integration.hpp
		include/kernels.hpp
		include/weak_magnetism.hpp
		include/opacities.hpp
		src/output/outputs.hpp
		include/eos.hpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/mnewt_GSL.cpp
		tests/tests_quadrature_integrate/gsl_routines/integration_GSL.cpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/kernel_nes.hpp
		tests/tests_opacities_m1/test_opacities_isoscatt.cpp
        tests/tests_opacities_m1/test_m1_opacities.cpp
		tests/tests.hpp
)
target_link_libraries(test_opacities_isoscatt Kokkos::kokkos)
target_link_libraries(test_opacities_isoscatt OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_opacities_isoscatt GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_opacities_isoscatt m)

# test opacities all reactions
add_executable(test_opacities_all_reactions
		include/constants.hpp
		include/bns_nurates.hpp
		include/distribution.hpp
		include/functions.hpp
		include/integration.hpp
		include/kernels.hpp
		include/weak_magnetism.hpp
		include/opacities.hpp
		src/output/outputs.hpp
		include/eos.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/kernel_nes.hpp
        tests/tests_opacities_m1/test_opacities_all_reactions.cpp
        tests/tests_opacities_m1/test_m1_opacities.cpp
		tests/tests.hpp
)
target_link_libraries(test_opacities_all_reactions Kokkos::kokkos)
target_link_libraries(test_opacities_all_reactions OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_opacities_all_reactions GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_opacities_all_reactions m)

add_executable(test_opacities_thc
		include/constants.hpp
		include/bns_nurates.hpp
		include/distribution.hpp
		include/functions.hpp
		include/integration.hpp
		include/kernels.hpp
		include/weak_magnetism.hpp
		include/opacities.hpp
		src/output/outputs.hpp
		include/eos.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/kernel_nes.hpp
		tests/tests_opacities_m1/test_opacities_thc.cpp
        tests/tests_opacities_m1/test_m1_opacities.cpp
		tests/tests.hpp
)
target_link_libraries(test_opacities_thc Kokkos::kokkos)
target_link_libraries(test_opacities_thc OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_opacities_thc GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_opacities_thc m)

# test m1 distribution
add_executable(test_m1_distribution
		include/constants.hpp
		include/bns_nurates.hpp
		include/distribution.hpp
		tests/tests_distribution/test_m1_distribution.cpp
)

target_link_libraries(test_m1_distribution Kokkos::kokkos)
target_link_libraries(test_m1_distribution OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_m1_distribution GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_m1_distribution m)

# test quadrature convergence spectral
add_executable(test_quadrature_convergence_spectral
		include/constants.hpp
		include/bns_nurates.hpp
		include/distribution.hpp
		include/functions.hpp
		include/integration.hpp
		include/kernels.hpp
		include/weak_magnetism.hpp
		include/opacities.hpp
		src/output/outputs.hpp
		include/eos.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/kernel_nes.hpp
		tests/tests_quadrature_convergence/test_quadrature_convergence_spectral.cpp
		tests/tests.hpp
)
target_link_libraries(test_quadrature_convergence_spectral Kokkos::kokkos)
target_link_libraries(test_quadrature_convergence_spectral OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_quadrature_convergence_spectral GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_quadrature_convergence_spectral m)

# test quadrature convergence integrated
add_executable(test_quadrature_convergence_integrated
		include/constants.hpp
		include/bns_nurates.hpp
		include/distribution.hpp
		include/functions.hpp
		include/integration.hpp
		include/kernels.hpp
		include/weak_magnetism.hpp
		include/opacities.hpp
		src/output/outputs.hpp
		include/eos.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/kernel_nes.hpp
		tests/tests_quadrature_convergence/test_quadrature_convergence_integrated.cpp
		tests/tests.hpp
)
target_link_libraries(test_quadrature_convergence_integrated Kokkos::kokkos)
target_link_libraries(test_quadrature_convergence_integrated OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_quadrature_convergence_integrated GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_quadrature_convergence_integrated m)

# test quadrature convergence find max
add_executable(test_quadrature_convergence_find_max
		include/constants.hpp
		include/bns_nurates.hpp
		include/distribution.hpp
		include/functions.hpp
		include/integration.hpp
		include/kernels.hpp
		include/weak_magnetism.hpp
		include/opacities.hpp
		src/output/outputs.hpp
		include/eos.hpp
		src/integration/integration_quadrature_schemes.cpp
		src/integration/integration_quadrature_io.cpp
		include/kernel_pair.hpp
		include/kernel_nes.hpp
		tests/tests_quadrature_convergence/test_quadrature_convergence_find_max.cpp
		tests/tests.hpp
)
target_link_libraries(test_quadrature_convergence_find_max Kokkos::kokkos)
target_link_libraries(test_quadrature_convergence_find_max OpenMP::OpenMP_CXX)
if (${GSL_FOUND})
	target_link_libraries(test_quadrature_convergence_find_max GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_quadrature_convergence_find_max m)


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
if (${GSL_FOUND})
	add_test(NAME test_pair_kernels
        	COMMAND $<TARGET_FILE:test_pair_kernels>)
	add_test(NAME test_pair_kernels_phi_mu_integrated
        	COMMAND $<TARGET_FILE:test_pair_kernels_phi_mu_integrated>)
endif()
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
