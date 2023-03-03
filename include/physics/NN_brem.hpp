#ifndef NN_BREM_H
#define NN_BREM_H

// \file NN_brem.hpp
// \brief Compute the nucleon-nucleon bremsstrahlung using the analytic fitting formula
//        in Hannestad & Raffelt 1998, Apj, 507, 339 (https://iopscience.iop.org/article/10.1086/306303/pdf)

#include <tuple>

#include "constants.hpp" // Header file for physical constants

using namespace std;
using namespace constants;

namespace brems {
  const double pi2   = pi*pi;
  const double pi1_2 = sqrt(pi);
  const double pi1_8 = pow(pi,1./8.);
  const double x_min   = 1.e-10;
  const double y_min   = 1.e-10;
  const double eta_min = 1.e-10;
}

/* Compute NN bremsstrahlung scattering kernel */
std::tuple<double,double> brem_kernel(const double w, const double wp, const double rho, const double tmev, const double xn, const double xp);

/* Compute s-component of the kernel */
double s_brem(double x, double y, double eta_star);

/* Compute g-component of the kernel */
double g_brem(double y, double eta_star);

/* Safe exp function to avoid overflow */
double fexp(const double x);

#endif
