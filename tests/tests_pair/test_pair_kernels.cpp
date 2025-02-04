//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_pair_kernels.c
//  \brief compare the PairKernels code with the Python implementation

#include <stdio.h>
#include <math.h>
#include <Kokkos_Core.hpp>

#include "kernel_pair.hpp"

using DevExeSpace  = Kokkos::DefaultExecutionSpace;
using DevMemSpace  = Kokkos::DefaultExecutionSpace::memory_space;
using HostMemSpace = Kokkos::HostSpace;
using LayoutWrapper =
    Kokkos::LayoutRight; // increments last index fastest views defined like

void test_pair_kernels()
{
    const double temp                       = 10.;
    const double mu_e                       = 84.181299999999993;
    const double pair_kernels_python[18][7] = {
        // nutype eta omega(MeV)  omega_prime(MeV) cos_theta R_p(cm^3s^-1)
        // R_a(cm^3s^-1)
        {0, 8.4181299999999997, 2, 2, 0, 3.7344281862494617e-37,
         5.5711121998146404e-37},
        {0, 8.4181299999999997, 2, 5, 0, 7.5952594997506983e-37,
         1.5294974381563826e-36},
        {0, 8.4181299999999997, 2, 50, 0, 7.6816308521376666e-37,
         1.3924664458243231e-34},
        {0, 8.4181299999999997, 5, 2, 0, 8.5842096421064313e-37,
         1.7286475408285996e-36},
        {0, 8.4181299999999997, 5, 5, 0, 1.7426501600901373e-36,
         4.7370142635342659e-36},
        {0, 8.4181299999999997, 5, 50, 0, 1.6703173370166646e-36,
         4.0871317668903464e-34},
        {0, 8.4181299999999997, 50, 2, 0, 3.6039543137528758e-36,
         6.5329687806960597e-34},
        {0, 8.4181299999999997, 50, 5, 0, 7.222632008133751e-36,
         1.7673197821036537e-33},
        {0, 8.4181299999999997, 50, 50, 0, 3.7405176147900196e-36,
         8.2390383297044382e-32},
        {1, 8.4181299999999997, 2, 2, 0, 7.9822789333374034e-38,
         1.1908160856214353e-37},
        {1, 8.4181299999999997, 2, 5, 0, 1.7098488384966924e-37,
         3.443212727887964e-37},
        {1, 8.4181299999999997, 2, 50, 0, 4.1187070226847197e-37,
         7.4660725562898844e-35},
        {1, 8.4181299999999997, 5, 2, 0, 1.7484861258538989e-37,
         3.5210186699128533e-37},
        {1, 8.4181299999999997, 5, 5, 0, 3.7248861049956034e-37,
         1.0125290212289139e-36},
        {1, 8.4181299999999997, 5, 50, 0, 8.4196498590303211e-37,
         2.0602203929942999e-34},
        {1, 8.4181299999999997, 50, 2, 0, 5.2266221327540312e-37,
         9.4744151143860744e-35},
        {1, 8.4181299999999997, 50, 5, 0, 1.0588883354811727e-36,
         2.5910143286093217e-34},
        {1, 8.4181299999999997, 50, 50, 0, 7.995294986861833e-37,
         1.7610809154750178e-32}};

    Kokkos::View<double*, LayoutWrapper, DevMemSpace> result_em_id_nue(
        "result_em_id_nue", 18);
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> result_abs_id_nue(
        "result_abs_id_nue", 18);
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> result_em_id_nux(
        "result_em_id_nux", 18);
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> result_abs_id_nux(
        "result_abs_id_nux", 18);

    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_result_em_id_nue(
        "h_result_em_id_nue", 18);
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_result_abs_id_nue(
        "h_result_abs_id_nue", 18);
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_result_em_id_nux(
        "h_result_em_id_nux", 18);
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_result_abs_id_nux(
        "h_result_abs_id_nux", 18);

    printf("=============================================================\n");
    printf("Pair process: Comparison of pair kernels with old Python code\n");
    printf("=============================================================\n");

    printf("omega omega' cos_theta R_em_e (new) R_em_e (python) %% err  "
           "R_abs_e (new)    R_abs_e (python) %% err  R_em_x (new)     R_em_x "
           "(python)  %% err  R_abs_x (new)    R_abs_x (python) %% err\n");

    Kokkos::parallel_for(
        "test_pair_kernels_loop", Kokkos::RangePolicy<>(DevExeSpace(), 0, 9),
        KOKKOS_LAMBDA(const int& i) {
            MyEOSParams eos_params              = {.temp = temp, .mu_e = mu_e};
            PairKernelParams pair_kernel_params = {
                .omega       = pair_kernels_python[i][2],
                .omega_prime = pair_kernels_python[i][3],
                .cos_theta   = pair_kernels_python[i][4],
                .mu          = -42.,
                .mu_prime    = -42.,
                .lmax        = 0,
                .filter      = 0};
            MyKernelOutput result =
                PairKernels(&eos_params, &pair_kernel_params);

            result_em_id_nue(i)  = result.em[id_nue];
            result_abs_id_nue(i) = result.abs[id_nue];
            result_em_id_nux(i)  = result.em[id_nux];
            result_abs_id_nux(i) = result.abs[id_nux];
        });

    Kokkos::deep_copy(h_result_em_id_nue, result_em_id_nue);
    Kokkos::deep_copy(h_result_abs_id_nue, result_abs_id_nue);
    Kokkos::deep_copy(h_result_em_id_nux, result_em_id_nux);
    Kokkos::deep_copy(h_result_abs_id_nux, result_abs_id_nux);

    for (int i = 0; i < 9; i++)
    {
        printf("%0.0e %0.0e %0.0e %0.10e %0.10e |%0.2f| %0.10e %0.10e |%0.2f| "
               "%0.10e %0.10e |%0.2f| %0.10e %0.10e |%0.2f|\n",
               pair_kernels_python[i][2], pair_kernels_python[i][3],
               pair_kernels_python[i][4], h_result_em_id_nue(i) * 1e-21,
               pair_kernels_python[i][5],
               Kokkos::fabs(
                   (h_result_em_id_nue(i) * 1e-21 - pair_kernels_python[i][5]) /
                   pair_kernels_python[i][5]) *
                   100.,
               h_result_abs_id_nue(i) * 1e-21, pair_kernels_python[i][6],
               Kokkos::fabs((h_result_abs_id_nue(i) * 1e-21 -
                             pair_kernels_python[i][6]) /
                                h_result_abs_id_nue(i) -
                            pair_kernels_python[i][6]) *
                   100,
               h_result_em_id_nux(i) * 1e-21, pair_kernels_python[i + 9][5],
               Kokkos::fabs((h_result_em_id_nux(i) * 1e-21 -
                             pair_kernels_python[i + 9][5]) /
                            pair_kernels_python[i + 9][5]) *
                   100.,
               h_result_abs_id_nux(i) * 1e-21, pair_kernels_python[i + 9][6],
               Kokkos::fabs((h_result_abs_id_nux(i) * 1e-21 -
                             pair_kernels_python[i + 9][6]) /
                            pair_kernels_python[i + 9][6]) *
                   100.);
    }
}

int main()
{
    Kokkos::initialize();
    test_pair_kernels();
    Kokkos::finalize();

    return 0;
}
