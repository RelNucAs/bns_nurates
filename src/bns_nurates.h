// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file bns_nurates.h
//  \brief essential data structures for the library

#ifndef BNS_NURATES_SRC_BNS_NURATES_H_
#define BNS_NURATES_SRC_BNS_NURATES_H_

// -------------------------------------------------
// Quadrature specific data structures

// Quadrature enum
// Holds the type of quadrature
enum Quadrature { kGauleg, kGaulag };
typedef enum Quadrature Quadrature;

// MyQuadrature struct
// Contains quadrature information, supports 1d/2d
struct MyQuadrature {
  enum Quadrature type;   // type of quadrature
  int dim;                // dimension of quadrature (1d/2d)
  int n;                  // number of points in the quadrature scheme
  double alpha;           // parameter for Gauss-Laguerre qaudrature (optional)
  double x1;              // lower limit of x
  double x2;              // upper limit of x
  double y1;              // lower limit of y (optional)
  double y2;              // upper limit of y (optional)
  double *x;              // points for the quadrature scheme (x)
  double *y;              // points for the quadrature scheme (y, optional)
  double *w;              // weights for the quadrature scheme
};
typedef struct MyQuadrature MyQuadrature;

// MyEOSParams struct
// Parameters which come from the EOS
// @TODO: decide if using standard or non-rel chemical potentials
struct MyEOSParams {
  double nb;      // baryon number density
  double temp;    // temperature
  double yp;      // proton fraction
  double yn;      // neutron fraction
  double mu_p;    // proton chemical potential
  double mu_n;    // neutron chemical potential
  double mu_e;    // electron chemical potential
  double mu_mu;   // muon chemical potential (this will be needed when including muon-dependent reactions)
  double dU;      // nuclear interaction correction on nuclear chemical potentials (as in Hempel 2015)
};
typedef struct MyEOSParams MyEOSParams;

// MyFunction struct
// A struct holding a function and its parameters
struct MyFunction {
  int dim;                                      // number of function variables (1/2)
  double (*function)(double var, void *params); // the function
  void *params;                                 // all parameters of the function
};
typedef struct MyFunction MyFunction;

// MyKernel struct
// Returns the absorption and production kernels for electron (e) and mu/tau (x) neutrinos
struct MyKernel {
  double production_e;
  double absorption_e;
  double production_x;
  double absorption_x;
  double production;
  double absorption;
};
typedef struct MyKernel MyKernel;

struct MyOpacityIntegrand {
  double em;
  double ab;
};
typedef struct MyOpacityIntegrand MyOpacityIntegrand;

// @FIXME: decide with Maitraya what to use for this
// Temporary struct for storing output of emissivity/absorptivity opacity output
// For the moment let's consider muonic (anti)neutrinos as nux
// Add in ab_(a)num, em_(a)num when we also muons will be considered
struct MyOpacity {
  double ab_nue;    // electron neutrino absorptivity
  double em_nue;    // electron neutrino emissivity
  double ab_anue;   // electron antineutrino absorptivity
  double em_anue;   // electron antineutrino emissivity
  //double ab_num;    // muon neutrino absorptivity
  //double em_num;    // muon neutrino emissivity
  //double ab_anum;   // muon antineutrino absorptivity
  //double em_anum;   // muon antineutrino emissivity
  double ab_nux;    // heavy (anti)neutrino absorptivity
  double em_nux;    // heavy (anti)neutrino emissivity
};
typedef struct MyOpacity MyOpacity;

#endif //BNS_NURATES_SRC_BNS_NURATES_H_
