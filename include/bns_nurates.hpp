// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file bns_nurates.h
//  \brief essential data structures for the library

#ifndef BNS_NURATES_SRC_BNS_NURATES_H_
#define BNS_NURATES_SRC_BNS_NURATES_H_

#include <stdbool.h>
#include <stddef.h>

// Define indices of neutrino species
#define id_nue 0
#define id_anue 1
#define id_nux 2
#define id_anux 3

#define total_num_species 4

// Define dimension of tabulated PairT function
#define dim_pair_t 10000

// Define maximum number of quadrature points
#define n_max 120

/* ==================================================================================
 * Integration structures
 * ==================================================================================
 */

/* Quadrature specific data structures
 *
 * Quadrature enum
 * Holds the type of quadrature: Gauss-Legendre, Gauss-Laguerre
 */
enum Quadrature
{
    kGauleg,
    kGaulag
};
typedef enum Quadrature Quadrature;

/* MyQuadrature struct
 *
 * Stores quadrature data and metadata, supports integration upto three
 * dimensions A default structure with metadata initialized is provided as
 * quadrature_default which specifies a 1d integration from 0 to 1 with
 * Gauss-Legendre and with 32 points. It is recommended that all MyQuadrature
 * data types are initialized with this before anything else is done with them:
 *
 * MyQuadrature quad = quadrature_default;
 *
 * Only the metadata is populated, so they have to be passed through the
 * necessary quadrature generation routines to populate the points and weights
 * arrays.
 *
 * Note: The weights and points array, irrespective of the number of dimensions
 * of the integration are always stored in a single 1d array. Any dimension
 * which is unused is populated with 1 in the weight and points arrays.
 *
 */
struct MyQuadrature
{
    enum Quadrature type; // type of quadrature (for the integration in the
                          // points variable, others are always kGauleg)
    double alpha;         // parameter for Gauss-Laguerre quadrature (optional)
    int dim;              // dimension of quadrature, can be 1,2,3
    int nx; // number of points in the quadrature scheme in the points direction
    int ny; // number of points in the quadrature scheme in the y direction, set
            // to 1 if not needed
    int nz; // number of points in the quadrature scheme in the z direction, set
            // to 1 if not needed
    double x1;      // lower limit of points, set to -42 if unused
    double x2;      // upper limit of points, set to -42 if unused
    double y1;      // lower limit of y, set to -42 if unused
    double y2;      // upper limit of y, set to -42 if unused
    double z1;      // lower limit of z, set to -42 if unused
    double z2;      // upper limit of z, set to -42 if unused
    double* points; // points for the quadrature scheme (store points in the
                    // points direction, then y and z in one flat array)
    double* w; // weights for the quadrature scheme (store points in the points
               // direction, then y and z in one flat array)
};
typedef struct MyQuadrature MyQuadrature;
__attribute__((unused)) static MyQuadrature quadrature_default = {.type =
                                                                      kGauleg,
                                                                  .alpha = 0.,
                                                                  .dim   = 1,
                                                                  .nx    = 20,
                                                                  .ny    = 1,
                                                                  .nz    = 1,
                                                                  .x1    = 0.,
                                                                  .x2    = 1.,
                                                                  .y1    = -42.,
                                                                  .y2    = -42.,
                                                                  .z1    = -42.,
                                                                  .z2 = -42.,
                                                                  .points = NULL,
                                                                  .w = NULL};

/* MyFunction struct
 *
 * A struct for holding one function and its parameters.
 * Use this when considering only one integrand.
 */
struct MyFunction
{
    int dim; // number of function variables (1/2)
    double (*function)(double* var, void* params); // function
    void* params;                                  // function parameters
};
typedef struct MyFunction MyFunction;

/* MyQuadratureIntegrand struct
 *
 * Holds metadata and integrand/integral data when multiple functions
 * are integrated in one go.
 *
 * This is used by the MyFunctionMultiD structure.
 */
struct MyQuadratureIntegrand
{
    int n;                // number of integrands/integrals
    double integrand[16]; // values of integrands/integrals (max: 16)
};
typedef struct MyQuadratureIntegrand MyQuadratureIntegrand;

/* MyFunctionMultiD struct
 *
 * A struct for holding multiple function and its parameters.
 * Use this when considering only multiple integrands.
 */
struct MyFunctionMultiD
{
    int dim; // number of function variables (1/2)
    MyQuadratureIntegrand (*function)(double* var3d,
                                      void* params); // the function
    MyQuadratureIntegrand
        my_quadrature_integrand; // integrand information and values
    void* params;                // function parameters
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
struct BremKernelParams
{
    double omega;            // neutrino energy before interaction
    double omega_prime;      // neutrino energy after interaction
    int l;                   // order of Legendre coefficient
    bool use_NN_medium_corr; // flag for inclusion of medium correction to HR98
                             // NN brem kernel as in Fischer16
};
typedef struct BremKernelParams BremKernelParams;

/* PairKernelParams struct
 *
 * Parameters for the pair kernel
 */
struct PairKernelParams
{
    double omega;       // neutrino energy
    double omega_prime; // anti-neutrino energy
    double cos_theta;   // cosine of angle between nu and a-nu
    double mu;          // cosine of neutrino polar angle
    double mu_prime;    // cosine of anti-neutrino polar angle
    double lmax;        // maximum value of l for Legendre expansion
    double filter;      // filter parameter for pair kernel positivity
    double alpha[dim_pair_t];
    double pair_t[6][dim_pair_t];
};
typedef struct PairKernelParams PairKernelParams;

/* PairKernelParams struct
 *
 * Parameters for the inelastic NES/NPS kernel
 */
struct InelasticScattKernelParams
{
    double omega;       // neutrino energy
    double omega_prime; // anti-neutrino energy
                        // @TODO: complete here
};
typedef struct InelasticScattKernelParams InelasticScattKernelParams;

/* MyKernelParams struct
 *
 * Unified structure for holding parameters for multiple kernels
 */
struct MyKernelParams
{
    PairKernelParams pair_kernel_params; // pair kernel parameters
    BremKernelParams brem_kernel_params; // Bremsstrahlung kernel
    InelasticScattKernelParams
        inelastic_kernel_params; // inelastic scattering kernel parameters
};
typedef struct MyKernelParams MyKernelParams;

// MyKernelOutput struct
//
/* MyKernelOutput struct
 *
 * holds emission/absorption related quantities for the pair and Bremsstrahlung
 * process and for inelastic scattering on leptons
 */
struct MyKernelOutput
{
    double em[total_num_species];  // emission/production kernel
    double abs[total_num_species]; // absorption/annihilation kernel
};
typedef struct MyKernelOutput MyKernelOutput;

// @TODO: decide what to do with the following
struct MyKernelQuantity
{
    double
        em_e; // quantity related to emission/production for electron neutrinos
    double abs_e; // quantity related to absorption for electron neutrinos
    double em_x; // quantity related to emission/production for mu/tau neutrinos
    double abs_x; // quantity related to absorption for mu/tau neutrinos
};
typedef struct MyKernelQuantity MyKernelQuantity;

/* ==================================================================================
 * EOS structures
 * ==================================================================================
 */

/* MyEOSParams struct
 *
 * Holds EOS parameters
 *
 * Chemical potentials include the rest mass contribution
 */
struct MyEOSParams
{
    double nb;    // baryon number density
    double temp;  // temperature
    double ye;    // electron fraction
    double yp;    // proton fraction
    double yn;    // neutron fraction
    double mu_p;  // proton chemical potential
    double mu_n;  // neutron chemical potential
    double mu_e;  // electron chemical potential
    double mu_mu; // muon chemical potential (this will be needed when including
                  // muon-dependent reactions)
    double dU; // nonrelativistic mean field intereaction potential difference
               // (as in Hempel 2015, Oertel et al. 2020)
    double
        dm_eff; // nonrelativistic mean field effective nucleon mass differenece
};
typedef struct MyEOSParams MyEOSParams;

/* ==================================================================================
 * Opacity structures
 * ==================================================================================
 */

/* MyOpacity struct
 *
 * Structure for storing emissivity and absorptivity output for
 * [0] electron neutrino (nue) [1] anti-electron neutrino (anue) [2] mu/tau
 * neutrino (nux)
 *
 * Can be used in the general case or the energy integrated case.
 *
 * Note: For the moment let's consider muonic (anti)neutrinos as nux
 *
 */
struct MyOpacity
{
    double abs[total_num_species]; // absorptivity
    double em[total_num_species];  // emissivity
};
typedef struct MyOpacity MyOpacity;

/* OpacityParams struct
 *
 * Store additional flags when computing opacities
 */
struct OpacityParams
{
    bool use_dU;     // flag for dU correction
    bool use_dm_eff; // flag for dm_eff correction
    bool use_WM_ab;  // flag for WM correction (and related) on absorption rates
    bool use_WM_sc;  // flag for WM correction (and related) on scattering rates
    bool use_decay;  // flag for inclusion of nucleon decay rates
    bool use_BRT_brem; // flag for computing NN brem rates using BRT06 instead
                       // of HR98
    bool use_NN_medium_corr; // flag for inclusion of medium correction to HR98
                             // NN brem kernel as in Fischer16
    bool neglect_blocking;   // flag for neglecting blocking factor of
                             // antineutrino in pair (nu + anu) processes
};
typedef struct OpacityParams OpacityParams;
__attribute__((unused)) static OpacityParams opacity_params_default_all = {
    .use_dU             = true,
    .use_dm_eff         = true,
    .use_WM_ab          = true,
    .use_WM_sc          = true,
    .use_decay          = true,
    .use_BRT_brem       = true,
    .use_NN_medium_corr = true,
    .neglect_blocking   = true};
__attribute__((unused)) static OpacityParams opacity_params_default_none = {
    .use_dU             = false,
    .use_dm_eff         = false,
    .use_WM_ab          = false,
    .use_WM_sc          = false,
    .use_decay          = false,
    .use_BRT_brem       = false,
    .use_NN_medium_corr = false,
    .neglect_blocking   = false};

/* ==================================================================================
 * M1 structures
 * ==================================================================================
 */

/* NuDistributionParams struct
 *
 * Structure for storing the parameters of the distribution function for
 * different species of neutrinos in the optically thin and thick regime.
 *
 * Supports [0]: electron neutrino, [1]: electron anti-neutrino, [2]: mu/tau
 * (anti)neutrino
 *
 * Follows the definition in Federico's notes
 */
struct NuDistributionParams
{

    // parameters for optically thick distribution function
    double w_t[total_num_species];    // contribution factor
    double temp_t[total_num_species]; // temperature
    double eta_t[total_num_species];  // degeneracy parameter

    // parameters for optically thin distribution function
    double w_f[total_num_species];    // contribution factor
    double temp_f[total_num_species]; // temperature
    double c_f[total_num_species];    // constant in power from Federico's notes
    double beta_f[total_num_species]; // from Federico's notes
};
typedef struct NuDistributionParams NuDistributionParams;

/* M1Quantities struct
 *
 * Stores the radiation number density, energy density and radiation flux for
 * all neutrino species from M1 Also stores the Eddington factor from the
 * closure
 */
struct M1Quantities
{
    double n[total_num_species];    // radiation number density
    double J[total_num_species];    // radiation energy density
    double H[total_num_species][4]; // radiation flux
    double chi[total_num_species];  // Eddington factor
};
typedef struct M1Quantities M1Quantities;

struct OpacityFlags
{
    int use_abs_em;
    int use_pair;
    int use_brem;
    int use_inelastic_scatt;
    int use_iso;
};
typedef struct OpacityFlags OpacityFlags;
__attribute__((unused)) static OpacityFlags opacity_flags_default_all = {
    .use_abs_em          = 1,
    .use_pair            = 1,
    .use_brem            = 1,
    .use_inelastic_scatt = 1,
    .use_iso             = 1};
__attribute__((unused)) static OpacityFlags opacity_flags_default_none = {
    .use_abs_em          = 0,
    .use_pair            = 0,
    .use_brem            = 0,
    .use_inelastic_scatt = 0,
    .use_iso             = 0};

/* GreyOpacityParams struct
 *
 * Stores all parameters needed for computing grey source coefficients for M1
 */
struct GreyOpacityParams
{
    OpacityParams opacity_pars; // spectral opacity input parameters
    MyKernelParams kernel_pars; // kernel input parameters
    MyEOSParams eos_pars;       // eos parameters
    NuDistributionParams
        distr_pars;             // neutrino distribution function parameters
    M1Quantities m1_pars;       // M1 related quantities
    OpacityFlags opacity_flags; // flags to turn on and off reactions
};
typedef struct GreyOpacityParams GreyOpacityParams;


/* M1Opacities struct
 *
 * Stores the emissivity, absorption and scattering coefficients
 * for electron neutrino (nue), electron anti-neutrino (anue) and mu/tau
 * neutrinos (nux) as in Radice et al. (2022)
 */
struct M1Opacities
{
    /* Number coefficients */
    double eta_0[total_num_species];     // number emissivity coefficient
    double kappa_0_a[total_num_species]; // number absorption coefficient

    /* Energy coefficients */
    double eta[total_num_species];     // energy emissivity coefficient
    double kappa_a[total_num_species]; // energy absorption coefficient
    double kappa_s[total_num_species]; // scattering coefficient
};
typedef struct M1Opacities M1Opacities;


/* M1Matrix struct
 *
 * Stores quantities related to the computation of M1 source
 * coefficients in 2D matrix form
 */

struct M1Matrix
{
    double** m1_mat_em[total_num_species];
    double** m1_mat_ab[total_num_species];
};
typedef struct M1Matrix M1Matrix;

struct M1MatrixKokkos2D
{
    double m1_mat_em[total_num_species][n_max][n_max];
    double m1_mat_ab[total_num_species][n_max][n_max];
};
typedef struct M1MatrixKokkos2D M1MatrixKokkos2D;

struct M1MatrixKokkos1D
{
    double m1_mat[12][n_max];
};
typedef struct M1MatrixKokkos1D M1MatrixKokkos1D;

/* SpectralOpacities struct
 *
 * Stores the emissivity and inverse mean free path
 * for electron neutrino (nue), electron anti-neutrino (anue) and mu/tau
 * neutrinos (nux)
 */
struct SpectralOpacities
{
    double j[total_num_species];       // emissivity
    double kappa[total_num_species];   // absorptivity / inverse mean free path
    double j_s[total_num_species];     // "emissivity" for scattering
    double kappa_s[total_num_species]; // inverse mean free path for scattering
};
typedef struct SpectralOpacities SpectralOpacities;

/* ==================================================================================
 * Legacy and unused structures @TODO: repurpose or remove safely
 * ==================================================================================
 */

// special function struct
// function returns 4 values
struct MyFunctionSpecial
{
    int dim; // number of function variables (1/2)
    MyKernelQuantity (*function)(
        double* var, MyEOSParams* my_eos_params,
        MyKernelParams* my_kernel_params); // the function
    MyEOSParams* eos_params;               // all eos parameters of the function
    MyKernelParams* kernel_params;         // all other parameters
};
typedef struct MyFunctionSpecial MyFunctionSpecial;

#endif // BNS_NURATES_SRC_BNS_NURATES_H_
