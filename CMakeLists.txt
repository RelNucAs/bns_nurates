cmake_minimum_required(VERSION 3.2)
project(bns_nurates VERSION 0.1 LANGUAGES C CXX)

# Kokkos require C++17 standard as minimum
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

# Force flag to denote C++ standard, as deprecated c++1z is not supported by nvcc.
set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++17")

# Kokkos recommends NOT using default compiler extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# following policy required to enable target_link_libraries() below
cmake_policy(SET CMP0079 NEW)

find_package(GSL)

# set OpenMP macro to OFF my default
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

set(BUILD_NURATES_THC OFF CACHE BOOL "Build nurates for THC")

include_directories(include)
include_directories(src)

################################################################################
add_executable(point_check
			include/constants.hpp
			include/bns_nurates.hpp
			include/distribution.hpp
			include/functions.hpp
			include/integration.hpp
			include/kernels.hpp
			include/weak_magnetism.hpp
			include/kernel_beta_iso.hpp
			src/output/outputs.hpp
			include/eos.hpp
			src/integration/integration_quadrature_io.cpp
			include/kernel_pair.hpp
			include/kernel_nes.hpp
      point_check.cpp
)
if(ENABLE_KOKKOS)
		target_compile_definitions(point_check PRIVATE KOKKOS_FLAG=1)
		target_link_libraries(point_check Kokkos::kokkos)
endif()
if (ENABLE_OPENMP)
		target_link_libraries(point_check OpenMP::OpenMP_CXX)
endif()
if (${GSL_FOUND})
		target_link_libraries(point_check GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(point_check m)
################################################################################


################################################################################
add_executable(test_opacities_benchmarks
			include/constants.hpp
			include/bns_nurates.hpp
			include/distribution.hpp
			include/functions.hpp
			include/integration.hpp
			include/kernels.hpp
			include/weak_magnetism.hpp
			include/kernel_beta_iso.hpp
			src/output/outputs.hpp
			include/eos.hpp
			src/integration/integration_quadrature_io.cpp
			include/kernel_pair.hpp
			include/kernel_nes.hpp
      tests/tests_opacities_m1/test_opacities_benchmarks.cpp
)
if(ENABLE_KOKKOS)
		target_compile_definitions(test_opacities_benchmarks PRIVATE KOKKOS_FLAG=1)
		target_link_libraries(test_opacities_benchmarks Kokkos::kokkos)
endif()
if (ENABLE_OPENMP)
		target_link_libraries(test_opacities_benchmarks OpenMP::OpenMP_CXX)
endif()
if (${GSL_FOUND})
		target_link_libraries(test_opacities_benchmarks GSL::gsl GSL::gslcblas)
endif()
target_link_libraries(test_opacities_benchmarks m)
################################################################################
