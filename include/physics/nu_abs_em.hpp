#ifndef NUABSEM_H
#define NUABSEM_H

// \file nu_abs_em.hpp
// \brief Computation of emissivity and absorptivity for neutrino absorption on neutron
//        (and inverse) and for antineutrino absorption on proton (and inverse)
//        Ref: Bruenn, 1985 (https://articles.adsabs.harvard.edu/pdf/1985ApJS...58..771B)

#include "constants.hpp" //Header file for physical constants

using namespace constants;

namespace nuabsem {
  const double g1 = (GF*GF/pi) / pow(h*c/(2.*pi),4.);
  const double g2 = gV*gV+3.*gA*gA;
  const double mu_thres = 1.e-2;
}

/* Factors resulting from nucleon phase space integration */
double eta_np(const double np, const double nn, const double mu_hat, const double T);
double eta_pn(const double np, const double nn, const double mu_hat, const double T);

/* Theta step function */
double theta(const double x);

/* Neutrino absorption on neutron (nul + n -> l- + p) */
std::tuple<double,double> nu_n_abs(const double omega, const double nb, const double T, const double lep_mass, const double yp, const double yn, const double mu_l, const double mu_hat, const double deltaU);
std::tuple<double,double> nu_n_abs_stim(const double omega, const double nb, const double T, const double lep_mass, const double yp, const double yn, const double mu_l, const double mu_hat, const double deltaU);

/* Antineutrino absorption on neutron (anul + p -> l+ + n) */
std::tuple<double,double> nu_p_abs(const double omega, const double nb, const double T, const double lep_mass, const double yp, const double yn, const double mu_l, const double mu_hat, const double deltaU);
std::tuple<double,double> nu_p_abs_stim(const double omega, const double nb, const double T, const double lep_mass, const double yp, const double yn, const double mu_l, const double mu_hat, const double deltaU);
