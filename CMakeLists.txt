cmake_minimum_required(VERSION 3.2)
project(bns_nurates VERSION 1.0 LANGUAGES C CXX)

# Kokkos require C++17 standard as minimum
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Og")

# Force flag to denote C++ standard, as deprecated c++1z is not supported by nvcc.
set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++17")

# Kokkos recommends NOT using default compiler extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# following policy required to enable target_link_libraries() below
cmake_policy(SET CMP0079 NEW)

find_package(GSL)

# set OpenMP macro to OFF
set(ENABLE_OPENMP ON CACHE BOOL "Enable OpenMP support")

if (ENABLE_OPENMP)
    message(STATUS "ENABLE_OPENMP set to ${ENABLE_OPENMP}: Compiling with OpenMP package.")
    find_package(OpenMP COMPONENTS CXX)
    if (OpenMP_CXX_FOUND)
        message(STATUS "OpenMP package found.")
    else()
        message(FATAL_ERROR "OpenMP package is required but could not be found.")
    endif()
else()
    message(STATUS "ENABLE_OPENMP set to ${ENABLE_OPENMP}: Not compiling with OpenMP package.")
endif()

set(ENABLE_KOKKOS ON CACHE BOOL "Enable Kokkos support")
if (ENABLE_KOKKOS)
    message(STATUS "ENABLE_KOKKOS set to ${ENABLE_KOKKOS}: Compiling with Kokkos package.")
else()
    message(STATUS "ENABLE_KOKKOS set to ${ENABLE_KOKKOS}: Not compiling with Kokkos package.")
endif()

if(ENABLE_KOKKOS AND ENABLE_OPENMP)
    set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "Compile Kokkos with OpenMP")
endif()

if(ENABLE_KOKKOS)
    add_subdirectory(kokkos)
    include_directories(${Kokkos_INCLUDE_DIRS_RET})
endif ()


if (ENABLE_OPENMP)
  set(OPENMP_PARALLEL_ENABLED 1)
else()
  set(OPENMP_PARALLEL_ENABLED 0)
endif()


if (${GSL_FOUND})
    add_compile_definitions(GSL_INCLUDES_H_)
endif()

add_definitions(-DSOURCE_DIR="${CMAKE_SOURCE_DIR}")

include_directories(include)
include_directories(tests)

################################################################################
macro(ADD_BNSNURATES_EXECUTABLE NAME PATH)
add_executable(
    ${NAME}
    ${PATH})
if(ENABLE_KOKKOS)
    target_compile_definitions(${NAME} PRIVATE KOKKOS_FLAG=1)
    target_link_libraries(${NAME} Kokkos::kokkos)
endif()
if (ENABLE_OPENMP)
    target_link_libraries(${NAME} OpenMP::OpenMP_CXX)
endif()
if (${GSL_FOUND})
    target_link_libraries(${NAME} GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(${NAME} m)
endmacro()
################################################################################

ADD_BNSNURATES_EXECUTABLE(
MWE
mwe.cpp)

ADD_BNSNURATES_EXECUTABLE(
point_out
utils/point_out.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_opacities_benchmarks
tests/tests_opacities_m1/test_opacities_benchmarks.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_beta
tests/tests_beta/test_beta.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_iso
tests/tests_iso/test_iso.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_brem
tests/tests_brem/test_brem.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_pair_kernels
tests/tests_pair/test_pair_kernels.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_gauss_legendre_1d
tests/tests_quadrature_integrate/test_gauss_legendre_1d.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_opacities_electron_capture
tests/tests_opacities_m1/test_opacities_electron_capture.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_opacities_isoscatt
tests/tests_opacities_m1/test_opacities_isoscatt.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_opacities_brem
tests/tests_opacities_m1/test_opacities_brem.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_opacities_pair
tests/tests_opacities_m1/test_opacities_pair.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_opacities_inelastic_scatt
tests/tests_opacities_m1/test_opacities_inelastic_scatt.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_opacities_all_reactions
tests/tests_opacities_m1/test_opacities_all_reactions.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_dumb_function
tests/tests_opacities_m1/test_dumb_function.cpp)
ADD_BNSNURATES_EXECUTABLE(
test_spectral_opacities_all_reactions
tests/tests_opacities_spectral/test_opacities_spectral.cpp)

ADD_BNSNURATES_EXECUTABLE(
test_m1_distribution
tests/tests_distribution/test_m1_distribution.cpp)
