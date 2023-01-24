//Opacities

#pragma once //compile only once

#include "constants.hpp"

using namespace constants;

namespace nuabsem
{

	const double g1 = (GF*GF/pi) / pow(h*c/(2.*pi),4.);
	const double g2 = gV*gV+3.*gA*gA;
}

double eta_np(const double np, const double nn, const double mu_hat, const double T);

double eta_pn(const double np, const double nn, const double mu_hat, const double T);

double theta(const double x);

//electron neutrino absorption on neutrons
std::tuple<double,double> nu_n_abs(const double omega, const double nb, const double T, const double lep_mass, const double yp, const double yn, const double mu_l, const double mu_hat, const double deltaU);

//electron antineutrino absorption on protons
std::tuple<double,double> nu_p_abs(const double omega, const double nb, const double T, const double lep_mass, const double yp, const double yn, const double mu_l, const double mu_hat, const double deltaU);
