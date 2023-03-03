#ifndef WEAK_MAG_H
#define WEAK_MAG_H

// !\file weak_magnetism.hpp
// \brief Evaluation of phase space/recoil/weak magnetism correction for (anti)neutrino
//        emission/absorption on nucleons and elastic scattering on nucleons
//        Ref: Horowitz, 2002 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)

/* Correction for electron neutrino absorption on neutron (nue + n -> e- + p) */
double WM_nue_abs(const double e_nu);

/* Correction for electron antineutrino absorption on proton (anue + p -> e+ + n) */
double WM_anue_abs(const double e_nu);

/* Correction for (anti)neutrino scattering on nucleons (nu + N -> nu + N) */
std::tuple<double,double> WM_scatt(const double Enu, const int reacflag);

#endif
