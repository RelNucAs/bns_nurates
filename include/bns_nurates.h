// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file bns_nurates.h
//  \brief essential data structures for the library

#ifndef BNS_NURATES_SRC_BNS_NURATES_H_
#define BNS_NURATES_SRC_BNS_NURATES_H_

#include <stdbool.h>

#define total_num_species 3

/* ==================================================================================
 * Integration structures
 * ==================================================================================
 */

/* Quadrature specific data structures
 *
 * Quadrature enum
 * Holds the type of quadrature: Gauss-Legendre, Gauss-Laguerre
 */
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
static MyQuadrature quadrature_default = {.type=kGauleg, .dim=1, .nx=55, .ny=1, .nz=1, .alpha=0., .x1=0., .x2=1., .y1=-42., .y2=-42., .z1=-42., .z2=-42.};

/* MyFunction struct
 *
 * A struct for holding one function and its parameters.
 * Use this when considering only one integrand.
 */
struct MyFunction {
  int dim;                                          // number of function variables (1/2)
  double (*function)(double *var, void *params);    // function
  void *params;                                     // function parameters
};
typedef struct MyFunction MyFunction;

/* MyQuadratureIntegrand struct
 *
 * Holds metadata and integrand/integral data when multiple functions
 * are integrated in one go.
 *
 * This is used by the MyFunctionMultiD structure.
 */
struct MyQuadratureIntegrand {
  int n;                  // number of integrands/integrals
  double integrand[12];   // values of integrands/integrals (max: 12)
};
typedef struct MyQuadratureIntegrand MyQuadratureIntegrand;

/* MyFunctionMultiD struct
 *
 * A struct for holding multiple function and its parameters.
 * Use this when considering only multiple integrands.
 */
struct MyFunctionMultiD {
  int dim;                                                        // number of function variables (1/2)
  MyQuadratureIntegrand (*function)(double *var3d, void *params); // the function
  MyQuadratureIntegrand my_quadrature_integrand;                  // integrand information and values
  void *params;                                                   // function parameters
};
typedef struct MyFunctionMultiD MyFunctionMultiD;

/* ==================================================================================
 * Kernel structures
 * ==================================================================================
 */

/* BremKernelParams struct
 *
 * Parameters for Bremsstrahlung kernel
 */
struct BremKernelParams {
  double omega;         // neutrino energy before interaction
  double omega_prime;   // neutrino energy after interaction
};
typedef struct BremKernelParams BremKernelParams;

/* PairKernelParams struct
 *
 * Parameters for the pair kernel
 */
struct PairKernelParams {
  double omega;         // neutrino energy
  double omega_prime;   // anti-neutrino energy
  double cos_theta;     // cosine of angle between nu and a-nu
  double mu;            // cosine of neutrino polar angle
  double mu_prime;      // cosine of anti-neutrino polar angle
  double lmax;          // maximum value of l for Legendre expansion
  double filter;        // filter parameter for pair kernel positivity
};
typedef struct PairKernelParams PairKernelParams;

/* MyKernelParams struct
 *
 * Unified structure for holding parameters for multiple kernels
 */
struct MyKernelParams {
  PairKernelParams pair_kernel_params;    // pair kernel parameters
  BremKernelParams brem_kernel_params;    // Bremsstrahlung kernel parameters
};
typedef struct MyKernelParams MyKernelParams;

// MyKernelQuantity struct
//
/* MyKernelQuantity struct
 *
 * holds emission/absorption related quantities for the pair and Bremsstrahlung process
 * for electron 'e' neutrinos and mu/tau 'x' neutrinos.
 */
struct MyKernelQuantity {
  double em_e;    // quantity related to emission/production for electron neutrinos
  double abs_e;   // quantity related to absorption for electron neutrinos
  double em_x;    // quantity related to emission/production for mu/tau neutrinos
  double abs_x;   // quantity related to absorption for mu/tau neutrinos
};
typedef struct MyKernelQuantity MyKernelQuantity;

/* ==================================================================================
 * EOS structures
 * ==================================================================================
 */

/* MyEOSParams struct
 *
 * Holds EOS parameters
 * @TODO: decide if using standard or non-rel chemical potentials
 */
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

/* ==================================================================================
 * Opacity structures
 * ==================================================================================
 */

/* MyOpacity struct
 *
 * Structure for storing emissivity and absorptivity output for
 * [0] electron neutrino (nue) [1] anti-electron neutrino (anue) [2] mu/tau neutrino (nux)
 *
 * Can be used in the general case or the energy integrated case.
 *
 * Note: For the moment let's consider muonic (anti)neutrinos as nux
 *       Add in ab_(a)num, em_(a)num when we also muons will be considered
 *
 * @TODO: Add in ab_(a)num, em_(a)num when muons are considered
 */
struct MyOpacity {
  double ab_nue;    // electron neutrino absorptivity
  double em_nue;    // electron neutrino emissivity
  double ab_anue;   // electron antineutrino absorptivity
  double em_anue;   // electron antineutrino emissivity
  double ab_nux;    // heavy (anti)neutrino absorptivity
  double em_nux;    // heavy (anti)neutrino emissivity
  //double ab_num;    // muon neutrino absorptivity
  //double em_num;    // muon neutrino emissivity
  //double ab_anum;   // muon antineutrino absorptivity
  //double em_anum;   // muon antineutrino emissivity
};
typedef struct MyOpacity MyOpacity;

/* OpacityParams struct
 *
 * Store additional flags when computing opacities
 */
struct OpacityParams {
  bool use_dU;     // flag for dU correction
  bool use_WM_ab;  // flag for WM correction (and related) on absorption rates
  bool use_WM_sc;  // flag for WM correction (and related) on scattering rates
};
typedef struct OpacityParams OpacityParams;

/* ==================================================================================
 * M1 structures
 * ==================================================================================
 */

/* NuDistributionParams struct
 *
 * Structure for storing the parameters of the distribution function for different species of neutrinos
 * in the optically thin and thick regime.
 *
 * Supports [0]: electron neutrino, [1]: electron anti-neutrino, [2]: mu/tau (anti)neutrino
 *
 * Follows the definition in Federico's notes
 */
struct NuDistributionParams {

  // parameters for optically thick distribution function
  double w_t[total_num_species];         // contribution factor
  double temp_t[total_num_species];      // temperature
  double eta_t[total_num_species];       // degeneracy parameter

  // parameters for optically thin distribution function
  double w_f[total_num_species];         // contribution factor
  double temp_f[total_num_species];      // temperature
  double c_f[total_num_species];         // constant in power from Ferederico's notes

};
typedef struct NuDistributionParams NuDistributionParams;

/* M1Quantities struct
 *
 * Stores the radiation number density, energy density and radiation flux for all neutrino species from M1
 * Also stores the Eddington factor from the closure
 */
struct M1Quantities {
  double n[total_num_species];      // radiation number density
  double J[total_num_species];      // radiation energy density
  double H[total_num_species][4];   // radiation flux
  double chi[total_num_species];    // Eddington factor
};
typedef struct M1Quantities M1Quantities;

struct OpacityFlags {
  int use_abs_em;
  int use_pair;
  int use_brem;
  int use_iso;
};
typedef struct OpacityFlags OpacityFlags;
static OpacityFlags opacity_flags_default_all = {.use_abs_em = 1., .use_pair = 1., .use_brem = 1., .use_iso = 1.};
static OpacityFlags opacity_flags_default_none = {.use_abs_em = 0., .use_pair = 0., .use_brem = 0., .use_iso = 0.};

/* GreyOpacityParams struct
 *
 * Stores all parameters needed for computing grey source coefficients for M1
 */
struct GreyOpacityParams {
  OpacityParams opacity_pars;      // spectral opacity input parameters
  MyKernelParams kernel_pars;      // kernel input parameters
  MyEOSParams eos_pars;            // eos parameters
  NuDistributionParams distr_pars; // neutrino distribution function parameters
  M1Quantities m1_pars;            // M1 related quantities
  OpacityFlags opacity_flags;      // flags to turn on and off reactions
};
typedef struct GreyOpacityParams GreyOpacityParams;

/* M1Opacities struct
 *
 * Stores the emissivity, absorption and scattering coefficients
 * for electron neutrino (nue), electron anti-neutrino (anue) and mu/tau neutrinos (nux)
 */
struct M1Opacities {
  /* Number coefficients */
  double eta_0_nue;       // number emissivity coefficient for nue
  double eta_0_anue;      // number emissivity coefficient for anue
  double eta_0_nux;       // number emissivity coefficient for nux

  double kappa_0_a_nue;   // number absorption coefficient for nue
  double kappa_0_a_anue;  // number absorption coefficient for anue
  double kappa_0_a_nux;   // number absorption coefficient for nux

  /* Energy coefficients */
  double eta_nue;       // energy emissivity coefficient for nue
  double eta_anue;      // energy emissivity coefficient for anue
  double eta_nux;       // energy emissivity coefficient for nux

  double kappa_a_nue;   // energy absorption coefficient for nue
  double kappa_a_anue;  // energy absorption coefficient for anue
  double kappa_a_nux;   // energy absorption coefficient for nux

  double kappa_s_nue;   // scattering coefficient for nue
  double kappa_s_anue;  // scattering coefficient for anue
  double kappa_s_nux;   // scattering coefficient for nux


};
typedef struct M1Opacities M1Opacities;

/* ==================================================================================
 * Legacy and unused structures @TODO: repurpose or remove safely
 * ==================================================================================
 */

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

// special function struct
// function returns 4 values
struct MyFunctionSpecial {
  int dim;                                                                                // number of function variables (1/2)
  MyKernelQuantity (*function)(double *var, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params);  // the function
  MyEOSParams *eos_params;                                                                // all eos parameters of the function
  MyKernelParams *kernel_params;                                                           // all other parameters
};
typedef struct MyFunctionSpecial MyFunctionSpecial;

#endif //BNS_NURATES_SRC_BNS_NURATES_H_
