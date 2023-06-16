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
        src/opacities/kernels/kernel_elastic_scatt.c
        src/opacities/weak_magnetism/weak_magnetism.h
        src/opacities/weak_magnetism/weak_magnetism.c
        src/opacities/weak_magnetism/nucfrmfac.c
	src/opacities/opacities.h
        src/opacities/nu_abs_em.c
        src/output/outputs.h
        src/eos/eos.h
        tests/tests.h
		tests/GSL_code/integration_GSL.c
		tests/GSL_code/mnewt_GSL.c
		tests/GSL_code/integration_GSL.c
		src/integration/integration_quadrature_schemes.c
		src/integration/integration_quadrature_io.c src/integration/integration_quadrature_integrate.c tests/tests_quadrature.c src/opacities/kernels/kernel_bremsstrahlung.c src/opacities/kernels/kernel_pair.c tests/tests_integration.c)

target_link_libraries(bns_nurates GSL::gsl GSL::gslcblas)
target_link_libraries(bns_nurates m)
