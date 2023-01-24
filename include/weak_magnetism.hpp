// Calculation of weak magnetism correction for neutrino emission/absorption and scattering

#pragma once //compile only once

double WM_nue_abs(const double e_nu);

double WM_anue_abs(const double e_nu);

std::tuple<double,double> WM_scatt(const double Enu, const int reacflag);
