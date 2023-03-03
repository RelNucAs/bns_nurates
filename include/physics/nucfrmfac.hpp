#ifndef NUCFRMFAC_H
#define NUCFRMFAC_H

//! \file nucfrmfac.hpp
//  \brief Calculation of single nucleon form factors as in C.J. Horowitz, 2002
//         (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001).
//         These are needed to compute the recoil and weak magnetism correction
//         for (anti)neutrino absorption on nucleons and elastic scattering on 
//         nucleons.

// Reactions are distinguished using the following indices:
//   reacflag = 0: (anti)neutrino scattering on proton  (nu p -> nu p)
//   reacflag = 1: (anti)neutrino scattering on neutron (nu n -> nu n)
//   reacflag = 2: (anti)neutrino absorption on nucleon (nue n -> e- p, anue p -> e+ n)

namespace formfactors {
  const double lamp =  1.793; //  proton magnetic moment?
  const double lamn = -1.913; // neutron magnetic moment
}

// Computation of single nucleon form factors for a given (anti)neutrino energy
// Inputs --> E: (anti)neutrino energy [MeV]
//            reacflag: index defining the reaction (see above)
std::tuple<double,double,double> nucfrmfac(const double E, const int reacflag);

#endif
