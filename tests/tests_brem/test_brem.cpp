//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_brem.c
//  \brief test routines for the bremsstrahlung process

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <Kokkos_Core.hpp>

#include "constants.hpp"
#include "bns_nurates.hpp"
#include "kernel_brem.hpp"

using DevExeSpace = Kokkos::DefaultExecutionSpace;
using DevMemSpace = Kokkos::DefaultExecutionSpace::memory_space;
using HostMemSpace = Kokkos::HostSpace;
using LayoutWrapper = Kokkos::LayoutRight;                // increments last index fastest views defined like

// Write data for plotting Fig. 3 of Hannestadt and Raffelt
void TestBremKernelS() {

  const int n = 100;
  double x[n];

  x[0] = 0.;
  x[n - 1] = 16.;
  double h = (x[n - 1] - x[0]) / ((double) n - 1);

  for (int i = 1; i < n - 1; i++) {
    x[i] = x[0] + i * h;
  }

  Kokkos::View<double*, LayoutWrapper, DevMemSpace> brem_kernel_s0("brem_kernel_s0", n);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> brem_kernel_s1("brem_kernel_s1", n);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> brem_kernel_s2("brem_kernel_s2", n);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> brem_kernel_s3("brem_kernel_s3", n);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> brem_kernel_s4("brem_kernel_s4", n);

  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_brem_kernel_s0("h_brem_kernel_s0", n);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_brem_kernel_s1("h_brem_kernel_s1", n);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_brem_kernel_s2("h_brem_kernel_s2", n);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_brem_kernel_s3("h_brem_kernel_s3", n);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_brem_kernel_s4("h_brem_kernel_s4", n);

  Kokkos::parallel_for("test_brem_kernels_s", Kokkos::RangePolicy<>(DevExeSpace(), 0, n),
  KOKKOS_LAMBDA(const int &i) {
    brem_kernel_s0(i) = BremKernelS(x[i], 0., 0.31);
    brem_kernel_s1(i) = BremKernelS(x[i], 0., 1.01);
    brem_kernel_s2(i) = BremKernelS(x[i], 0., 2.42);
    brem_kernel_s3(i) = BremKernelS(x[i], 2., 0.31);
    brem_kernel_s4(i) = BremKernelS(x[i], 4., 0.31);
  });

  Kokkos::deep_copy(h_brem_kernel_s0, brem_kernel_s0);
  Kokkos::deep_copy(h_brem_kernel_s1, brem_kernel_s1);
  Kokkos::deep_copy(h_brem_kernel_s2, brem_kernel_s2);
  Kokkos::deep_copy(h_brem_kernel_s3, brem_kernel_s3);
  Kokkos::deep_copy(h_brem_kernel_s4, brem_kernel_s4);

  for (int i = 0; i < n; i++) {
    printf("%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n", x[i], h_brem_kernel_s0(i), h_brem_kernel_s1(i), h_brem_kernel_s2(i),
           h_brem_kernel_s3(i), h_brem_kernel_s4(i));
  }
}

// Write data for plotting Fig. 4 of Hannestadt and Raffelt
void TestBremKernelG() {

  const int n = 100;
  double x[n];

  x[0] = 0.;
  x[n - 1] = 8.;
  double h = (x[n - 1] - x[0]) / ((double) n - 1);

  for (int i = 1; i < n - 1; i++) {
    x[i] = x[0] + i * h;
  }

  Kokkos::View<double*, LayoutWrapper, DevMemSpace> brem_kernel_g0("brem_kernel_g0", n);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> brem_kernel_g1("brem_kernel_g1", n);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> brem_kernel_g2("brem_kernel_g2", n);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> brem_kernel_g3("brem_kernel_g3", n);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> brem_kernel_g4("brem_kernel_g4", n);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> brem_kernel_g5("brem_kernel_g5", n);

  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_brem_kernel_g0("h_brem_kernel_g0", n);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_brem_kernel_g1("h_brem_kernel_g1", n);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_brem_kernel_g2("h_brem_kernel_g2", n);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_brem_kernel_g3("h_brem_kernel_g3", n);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_brem_kernel_g4("h_brem_kernel_g4", n);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_brem_kernel_g5("h_brem_kernel_g5", n);

  Kokkos::parallel_for("test_brem_kernels_g", Kokkos::RangePolicy<>(DevExeSpace(), 0, n),
  KOKKOS_LAMBDA(const int &i) {
    brem_kernel_g0(i) = BremKernelG(x[i], 0.31);
    brem_kernel_g1(i) = BremKernelG(x[i], 1.01);
    brem_kernel_g2(i) = BremKernelG(x[i], 2.42);
    brem_kernel_g3(i) = BremKernelG(x[i], 4.22);
    brem_kernel_g4(i) = BremKernelG(x[i], 6.14);
    brem_kernel_g5(i) = BremKernelG(x[i], 8.11);
  });

  Kokkos::deep_copy(h_brem_kernel_g0, brem_kernel_g0);
  Kokkos::deep_copy(h_brem_kernel_g1, brem_kernel_g1);
  Kokkos::deep_copy(h_brem_kernel_g2, brem_kernel_g2);
  Kokkos::deep_copy(h_brem_kernel_g3, brem_kernel_g3);
  Kokkos::deep_copy(h_brem_kernel_g4, brem_kernel_g4);
  Kokkos::deep_copy(h_brem_kernel_g5, brem_kernel_g5);

  for (int i = 0; i < n; i++) {
    printf("%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n", x[i], h_brem_kernel_g0(i), h_brem_kernel_g1(i), h_brem_kernel_g2(i),
           h_brem_kernel_g3(i), h_brem_kernel_g4(i), h_brem_kernel_g5(i));
  }
}

int main() {

  Kokkos::initialize();

  printf("Test Brem Kernels S:\n");
  printf("\n");
  TestBremKernelS();
  Kokkos::fence();
  printf("\n");
  printf("Test Brem Kernels G:\n");
  TestBremKernelG();
  Kokkos::fence();
  printf("\n");

  const double ref_em_ker  = 6.29472e-30; // [cm^3 s^-1]
  const double ref_abs_ker = 7.60512e-29; // [cm^3 s^-1]

  printf("Test Brem Main:\n");
  Kokkos::parallel_for("test_brem_kernels", Kokkos::RangePolicy<>(DevExeSpace(), 0, 1),
                       KOKKOS_LAMBDA(const int &i) {
                         BremKernelParams kernel_pars;

                         kernel_pars.omega = 10.; // [MeV]
                         kernel_pars.omega_prime = 20.; // [MeV]
                         kernel_pars.l = 0;

                         // EOS parameters
                         MyEOSParams eos_pars;

                         // @TODO: decide a global baryon mass to convert number to mass density and vice versa
                         eos_pars.nb = 3.73e+14 / kMb; // [cm-3]
                         eos_pars.temp = 12.04;  // [MeV]
                         eos_pars.yn = 0.6866;
                         eos_pars.yp = 0.3134;

                         // Compute kernels
                         MyKernelOutput
                             out = BremKernelsLegCoeff(&kernel_pars, &eos_pars);

                         // Print kernels
                         printf("Output emission   kernel = %.5e cm^3/s\n",
                                out.em[id_nue]);
                         printf("Output absorption kernel = %.5e cm^3/s\n",
                                out.abs[id_nue]);
                         printf("\n");

                         // Compare kernels
                         const double tol = 1.0E-06; // tolerance on relative difference
                         const double em_diff = fabs(out.em[id_nue] - ref_em_ker)
                             / ref_em_ker;  // emission kernel difference
                         const double abs_diff = fabs(out.abs[id_nue] - ref_abs_ker)
                             / ref_abs_ker; // absorption kernel difference

                         if ((em_diff > tol) || (abs_diff > tol)) {
                           printf("Test failed!\n");
                           printf("Comparison values are:\n");
                           printf("Emission   kernel = %.5e cm^3/s\n", ref_em_ker);
                           printf("Absorption kernel = %.5e cm^3/s\n", ref_abs_ker);
                         } else {
                           printf("Test passed!\n");
                         }
                       });
  Kokkos::finalize();

  return 0;
}
