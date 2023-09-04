// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file bns_nurates.h
//  \brief essential data structures for the library

#ifndef BNS_NURATES_SRC_BNS_NURATES_H_
#define BNS_NURATES_SRC_BNS_NURATES_H_

#include <stdbool.h>

//#define kNSp 3 // number of neutrino species

// -------------------------------------------------
// Quadrature specific data structures

// Quadrature enum
// Holds the type of quadrature
enum Quadrature { kGauleg, kGaulag };
typedef enum Quadrature Quadrature;

/* MyQuadrature struct
 *
 * Stores quadrature data and metadata, supports integration upto three dimensions
 * A default structure with metadata initialized is provided as quadrature_default which specifies
 * a 1d integration from 0 to 1 with Gauss-Legendre and with 32 points. It is recommended that all
 * MyQuadrature data types are initialized with this before anything else is done with them:
 *
 * MyQuadrature quad = quadrature_default;
 *
 * Only the metadata is populated, so they have to be passed through the necessary quadrature generation
 * routines to populate the points and weights arrays.
 *
 * Note: The weights and points array, irrespective of the number of dimensions of the integration are always
 * stored in a single 1d array. Any dimension which is unused is populated with 1 in the weight and points arrays.
 *
 */
struct MyQuadrature {
  enum Quadrature type;   // type of quadrature (for the integration in the points variable, others are always kGauleg)
  double alpha;           // parameter for Gauss-Laguerre quadrature (optional)
  int dim;                // dimension of quadrature, can be 1,2,3
  int nx;                 // number of points in the quadrature scheme in the points direction
  int ny;                 // number of points in the quadrature scheme in the y direction, set to 1 if not needed
  int nz;                 // number of points in the quadrature scheme in the z direction, set to 1 if not needed
  double x1;              // lower limit of points, set to -42 if unused
  double x2;              // upper limit of points, set to -42 if unused
  double y1;              // lower limit of y, set to -42 if unused
  double y2;              // upper limit of y, set to -42 if unused
  double z1;              // lower limit of z, set to -42 if unused
  double z2;              // upper limit of z, set to -42 if unused
  double *points;         // points for the quadrature scheme (store points in the points direction, then y and z in one flat array)
  double *w;              // weights for the quadrature scheme (store points in the points direction, then y and z in one flat array)
};
typedef struct MyQuadrature MyQuadrature;
static MyQuadrature quadrature_default = {.type=kGauleg, .dim=1, .nx=32, .ny=1, .nz=1, .alpha=0., .x1=0., .x2=1., .y1=-42., .y2=-42., .z1=-42., .z2=-42.};

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

struct var3d {
  double x;
  double y;
  double z;
};
typedef struct var3d var3d;
static var3d var3d_default = {.x = 0., .y = 0., .z = 0.};

// MyFunction struct
// A struct holding a function and its parameters
struct MyFunction {
  int dim;                                      // number of function variables (1/2)
  double (*function)(double *var3d, void *params); // the function
  void *params;                                 // all parameters of the function
};
typedef struct MyFunction MyFunction;

// bremsstrahlung kernel specific parameters
struct BremKernelParams {
  double omega;
  double omega_prime;
};
typedef struct BremKernelParams BremKernelParams;

// pair kernel specific parameters
struct PairKernelParams {
  double omega;
  double omega_prime;
  double cos_theta;
  double mu;
  double mu_prime;
  double lmax;
  double filter;
};
typedef struct PairKernelParams PairKernelParams;

// unified kernel params
struct MyKernelParams {
  PairKernelParams pair_kernel_params;
  BremKernelParams brem_kernel_params;
};
typedef struct MyKernelParams MyKernelParams;

// MyKernel struct
// Returns the absorption and production kernels for electron (e) and mu/tau (points) neutrinos
struct MyKernel {
  double production_e;
  double absorption_e;
  double production_x;
  double absorption_x;
  double production;
  double absorption;
};
typedef struct MyKernel MyKernel;

// MyOpacityQuantity struct
// holds emission/absorption quantities for two neutrino species
struct MyOpacityQuantity {
  double em_e;    // quantity related to emission/production for e neutrinos
  double abs_e;   // quantity related to absorption for e neutrinos
  double em_x;    // quantity related to emission/production for points neutrinos
  double abs_x;   // quantity related to absorption for points neutrinos
};
typedef struct MyOpacityQuantity MyOpacityQuantity;

// special function struct
// function returns 4 values
struct MyFunctionSpecial {
  int dim;                                                                                // number of function variables (1/2)
  MyOpacityQuantity (*function)(double *var, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params);  // the function
  MyEOSParams *eos_params;                                                                // all eos parameters of the function
  MyKernelParams *kernel_params;                                                           // all other parameters
};
typedef struct MyFunctionSpecial MyFunctionSpecial;

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

// Opacity parameters
struct OpacityParams {
  bool use_dU;     // flag for dU correction
  bool use_WM_ab;  // flag for WM correction (and related) on absorption rates
  bool use_WM_sc;  // flag for WM correction (and related) on scattering rates
};
typedef struct OpacityParams OpacityParams;

// @TODO: generalize the following structure to all neutrino species (e.g. double w_t -> double w_t[Nsp])

// NuDistributionParams struct
// structure needed for storing the parameters of the 
// neutrino distribution function
struct NuDistributionParams {
  // optically thick
  double w_t;
  double temp_t;
  double eta_t;
  // optically thin
  double w_f;
  double temp_f;
  double c_f;
};
typedef struct NuDistributionParams NuDistributionParams;

// @TODO: generalize the following structure to all neutrino species (e.g. double n -> double n[Nsp])
struct M1Quantities {
  double n;    // radiation number density
  double J;    // radiation energy density
  double H[4]; // radiation flux components (only three are independent since
  //                            H^alpha u_alpha = 0)
  double chi;  // closure
};
typedef struct M1Quantities M1Quantities;

// GreyOpacityParams struct
// structure for storing the parameters needed for the 
// computation of grey source coefficients
struct GreyOpacityParams {
  OpacityParams opacity_pars;      // spectral opacity input parameters
  MyKernelParams kernel_pars;      // kernel input parameters
  MyEOSParams eos_pars;            // eos parameters
  NuDistributionParams distr_pars; // neutrino distribution function parameters
  M1Quantities m1_pars;            // M1 related quantities
};
typedef struct GreyOpacityParams GreyOpacityParams;

// SourceCoeffs struct
// structure needed for storing the values of the 
// emissivity and absoprtion/scattering coefficients
// for the different neutrino species as in Radice
// et al. (2022)
struct SourceCoeffs {
  double R_nue;    // number, electron-type neutrino
  double R_anue;   // number, electron-type antineutrino
  //double R_num;    // number, muon-type neutrino
  //double R_anum;   // number, muon-type antineutrino
  double R_nux;    // number, heavy-type (anti)neutrino
  double Q_nue;    // energy, electron-type neutrino
  double Q_anue;   // energy, electron-type antineutrino
  //double Q_num;    // energy, muon-type neutrino
  //double Q_anum;   // energy, muon-type antineutrino
  double Q_nux;    // energy, heavy-type (anti)neutrino
};
typedef struct SourceCoeffs SourceCoeffs;

struct M1Opacities {
  double eta_e;
  double kappa_a_e;
  double kappa_s_e;

  double eta_x;
  double kappa_a_x;
  double kappa_s_x;
};
typedef struct M1Opacities M1Opacities;

struct MyQuadratureIntegrand {
  int n;    // number of integrands (maximum is 6)
  double integrand[6];
};
typedef struct MyQuadratureIntegrand MyQuadratureIntegrand;

struct MyFunctionMultiD {
  int dim;                                                        // number of function variables (1/2)
  MyQuadratureIntegrand (*function)(double *var3d, void *params); // the function
  MyQuadratureIntegrand my_quadrature_integrand;
  void *params;                                                   // all parameters of the function
};
typedef struct MyFunctionMultiD MyFunctionMultiD;

#endif //BNS_NURATES_SRC_BNS_NURATES_H_
