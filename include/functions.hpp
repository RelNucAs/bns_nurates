//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_SRC_FUNCTIONS_FUNCTIONS_H_
#define BNS_NURATES_SRC_FUNCTIONS_FUNCTIONS_H_

#include "bns_nurates.hpp"
#include "constants.hpp"

constexpr BS_REAL zero                  = 0;
constexpr BS_REAL one_twentyth          = 0.05;
constexpr BS_REAL one_tenth             = 0.1;
constexpr BS_REAL one_fifth             = 0.2;
constexpr BS_REAL one_third             = 0.3333333333333333333;
constexpr BS_REAL one                   = 1;
constexpr BS_REAL two                   = 2;
constexpr BS_REAL five                  = 5;
constexpr BS_REAL ten                   = 10;
constexpr BS_REAL twenty                = 20;
constexpr BS_REAL forty                 = 40;
constexpr BS_REAL onethousandsixhundred = 1600;
constexpr BS_REAL fdi_litconst          = 7.38905609893065023;

// Exception handling from Numerical Recipes
// @TODO: decide how to handle errors in the code

/*
#ifndef _USENRERRORCLASS_
#define throw(message)                                                         \
    //} //\
        //BS_PRINTF("ERROR: %s\n     in file %s at line %d\n", message,
__FILE__, \
               __LINE__);
                                                                     \
        //exit(1);
        //BS_REAL dummy_var = -42.;
        BS_PRINTF("Throw function here");
                                                                      \
    //}
#else
struct NRerror
{
    char* message;
    char* file;
    int line;
    NRerror(char* m, char* f, int l) : message(m), file(f), line(l)
    {
    }
};
#define throw(message) throw(NRerror(message, __FILE__, __LINE__));
void NRcatch(NRerror err)
{
    BS_PRINTF("ERROR: %s\n     in file %s at line %d\n", err.message, err.file,
           err.line);
    exit(1);
}
#endif
*/

/*===========================================================================*/

// safe_exp.c

// Safe exp function to avoid underflow/overflow
KOKKOS_INLINE_FUNCTION
BS_REAL SafeExp(const BS_REAL x)
{
    return exp(fmin(fmax(x, kBS_ExpLowLim), kBS_ExpUppLim));
}

/*===========================================================================*/

// bessel.c

/*
Evaluate modified Bessel functions of order v = 1
Description:
              The differential equation

                       2
                   2  d w       dw      2   2
                  x . --- + x . --- - (x + v ).w = 0
                        2       dx
                      dx

              has two solutions called modified Bessel functions
              Iv(x) and Kv(x).
              The routines bessi1 and bessk1 return the I and K for
              v = 1.
*/

// static
KOKKOS_INLINE_FUNCTION
BS_REAL bessi1(const BS_REAL x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1.  */
/*------------------------------------------------------------*/
{
    BS_REAL ax, ans;
    BS_REAL y;

    constexpr BS_REAL zero            = 0;
    constexpr BS_REAL fifteen_fourths = 3.75;
    constexpr BS_REAL A[7] = {0.5,          0.87890594,  0.51498869, 0.15084934,
                              0.2658733e-1, 0.301532e-2, 0.32411e-3};
    constexpr BS_REAL B[9] = {0.2282967e-1, -0.2895312e-1, 0.1787654e-1,
                              0.420059e-2,  0.39894228,    -0.3988024e-1,
                              -0.362018e-2, 0.163801e-2,   -0.1031555e-1};


    if ((ax = fabs(x)) < fifteen_fourths)
    {
        y = x / fifteen_fourths, y = y * y;
        ans =
            ax * (A[0] +
                  y * (A[1] +
                       y * (A[2] +
                            y * (A[3] + y * (A[4] + y * (A[5] + y * A[6]))))));
    }
    else
    {
        y   = fifteen_fourths / ax;
        ans = B[0] + y * (B[1] + y * (B[2] - y * B[3]));
        ans =
            B[4] + y * (B[5] + y * (B[6] + y * (B[7] + y * (B[8] + y * ans))));
        ans *= (exp(ax) / sqrt(ax));
    }

    return x < zero ? -ans : ans;
}


// static
KOKKOS_INLINE_FUNCTION
BS_REAL bessk1(const BS_REAL x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=1.  */
/*------------------------------------------------------------*/
{
    BS_REAL y, ans;

    constexpr BS_REAL two  = 2;
    constexpr BS_REAL four = 4;
    constexpr BS_REAL A[8] = {1.0,          1.0,         0.15443144,
                              -0.67278579,  -0.18156897, -0.1919402e-1,
                              -0.110404e-2, -0.4686e-4};
    constexpr BS_REAL B[7] = {1.25331414,   0.23498619,   -0.3655620e-1,
                              0.1504268e-1, -0.780353e-2, 0.325614e-2,
                              -0.68245e-3};


    if (x <= two)
    {
        y = x * x / four;
        ans =
            (log(x / two) * bessi1(x)) +
            (A[0] / x) *
                (A[1] +
                 y * (A[2] +
                      y * (A[3] +
                           y * (A[4] + y * (A[5] + y * (A[6] + y * (A[7])))))));
    }
    else
    {
        y   = two / x;
        ans = (exp(-x) / sqrt(x)) *
              (B[0] +
               y * (B[1] +
                    y * (B[2] +
                         y * (B[3] + y * (B[4] + y * (B[5] + y * (B[6])))))));
    }

    return ans;
}

/* Numerical recipes code

const BS_REAL k1pi[]={0.5,5.598072040178741e-2,1.818666382168295e-3,
2.397509908859959e-5,1.239567816344855e-7};
const BS_REAL k1qi[]={9.870202601341150e-1,1.292092053534579e-2,
5.881933053917096e-5};
const BS_REAL k1p[]={-3.079657578292062e-1,-8.109417631822442e-2,
-3.477550948593604e-3,-5.385594871975406e-5,-3.110372465429008e-7};
const BS_REAL k1q[]={9.861813171751389e-1,1.375094061153160e-2,
6.774221332947002e-5};
const BS_REAL k1pp[]={1.253314137315502,1.457171340220454e1,
6.063161173098803e1,1.147386690867892e2,1.040442011439181e2,
4.356596656837691e1,7.265230396353690,3.144418558991021e-1};
const BS_REAL k1qq[]={1.0,1.125154514806458e1,4.427488496597630e1,
7.616113213117645e1,5.863377227890893e1,1.850303673841586e1,
1.857244676566022,2.538540887654872e-2};

BS_REAL k1(const BS_REAL x) {
    // Returns the modiﬁed Bessel function K1(x) for positive real x.
    if (x <= 1.0) {  // Use two rational approximations.
        const BS_REAL z=x*x;
        const BS_REAL term = poly(k1pi,4,z)*log(x)/poly(k1qi,2,1.-z);
        return x*(poly(k1p,4,z)/poly(k1q,2,1.-z)+term)+1./x;
    } else {         // Rational approximation with e^{-x}/sqrt(x) factored out.
        const BS_REAL z=1.0/x;
        return exp(-x)*poly(k1pp,7,z)/(poly(k1qq,7,z)*sqrt(x));
    }
}

inline BS_REAL poly(const BS_REAL *cof, const int n, const BS_REAL x) {
    // Common code: Evaluate a polynomial.
    BS_REAL ans = cof[n];
    for (int i=n-1;i>=0;i--) ans = ans*x+cof[i];
    return ans;
}
*/

/*===========================================================================*/

// digamma.c

/*
 * Computation of Psi (Digamma) function
 *  Adapted from GNU Scientific Library (GSL 2.7.1)
 */

/* Chebyshev fits from SLATEC code for psi(x)

 Series for PSI        on the interval  0.         to  1.00000D+00
                                       with weighted error   2.03E-17
                                        log weighted error  16.69
                              significant figures required  16.39
                                   decimal places required  17.37

 Series for APSI       on the interval  0.         to  2.50000D-01
                                       with weighted error   5.54E-17
                                        log weighted error  16.26
                              significant figures required  14.42
                                   decimal places required  16.86

*/

typedef struct
{
    BS_REAL val;
    BS_REAL err;
} SFResult;

struct cheb_series_struct
{
    BS_REAL* c;   /* coefficients                */
    int order;    /* order of expansion          */
    BS_REAL a;    /* lower interval point        */
    BS_REAL b;    /* upper interval point        */
    int order_sp; /* effective single precision order */
};

typedef struct cheb_series_struct ChebSeries;

static BS_REAL psics_data[23] = {
    -.038057080835217922, .491415393029387130,  -.056815747821244730,
    .008357821225914313,  -.001333232857994342, .000220313287069308,
    -.000037040238178456, .000006283793654854,  -.000001071263908506,
    .000000183128394654,  -.000000031353509361, .000000005372808776,
    -.000000000921168141, .000000000157981265,  -.000000000027098646,
    .000000000004648722,  -.000000000000797527, .000000000000136827,
    -.000000000000023475, .000000000000004027,  -.000000000000000691,
    .000000000000000118,  -.000000000000000020};

static ChebSeries psi_cs = {psics_data, 22, -1, 1, 17};

static BS_REAL apsics_data[16] = {
    -.0204749044678185, -.0101801271534859, .0000559718725387,
    -.0000012917176570, .0000000572858606,  -.0000000038213539,
    .0000000003397434,  -.0000000000374838, .0000000000048990,
    -.0000000000007344, .0000000000001233,  -.0000000000000228,
    .0000000000000045,  -.0000000000000009, .0000000000000002,
    -.0000000000000000};

static ChebSeries apsi_cs = {apsics_data, 15, -1, 1, 9};

// Evaluation of the Chebyshev series cs at a given point x
KOKKOS_INLINE_FUNCTION
void ChebEvalE(const ChebSeries* cs, const BS_REAL x, SFResult* result)
{
    int j;
    BS_REAL d  = 0.0;
    BS_REAL dd = 0.0;

    constexpr BS_REAL half            = 0.5;
    constexpr BS_REAL two             = 2;
    constexpr BS_REAL kBS_dbl_epsilon = DBL_EPSILON;

    BS_REAL y  = (two * x - cs->a - cs->b) / (cs->b - cs->a);
    BS_REAL y2 = two * y;

    BS_REAL e = 0.0;

    for (j = cs->order; j >= 1; j--)
    {
        BS_REAL temp = d;
        d            = y2 * d - dd + cs->c[j];
        e += fabs(y2 * temp) + fabs(dd) + fabs(cs->c[j]);
        dd = temp;
    }

    {
        BS_REAL temp = d;
        d            = y * d - dd + half * cs->c[0];
        e += fabs(y * temp) + fabs(dd) + half * fabs(cs->c[0]);
    }

    result->val = d;
    result->err = kBS_dbl_epsilon * e + fabs(cs->c[cs->order]);

    return; // GSL_SUCCESS;
}


// Evaluation of Psi (Digamma) function (result && error)
/* digamma for x both positive and negative; we do both
 * cases here because of the way we use even/odd parts
 * of the function
 */
/*
KOKKOS_INLINE_FUNCTION
void SFPsiOutput(const BS_REAL x, SFResult* result)
{
    const BS_REAL y = fabs(x);

    if (x == 0.0 || x == -1.0 || x == -2.0)
    {
        // result->val = std::nan;
        // result->err = std::nan;
        BS_PRINTF("Error in digamma computation\n");
        exit(EXIT_FAILURE);
    }
    else if (y >= 2.0)
    {
        const BS_REAL t = 8.0 / (y * y) - 1.0;
        SFResult result_c;
        ChebEvalE(&apsi_cs, t, &result_c);
        if (x < 0.0)
        {
            const BS_REAL s = sin(kBS_Pi * x);
            const BS_REAL c = cos(kBS_Pi * x);
            if (fabs(s) < 2.0 * sqrt(DBL_MIN))
            {
                // result->val = std::nan;
                // result->err = std::nan;
                BS_PRINTF("Error in digamma computation\n");
                exit(EXIT_FAILURE);
            }
            else
            {
                result->val = log(y) - 0.5 / x + result_c.val - kBS_Pi * c / s;
                result->err = kBS_Pi * fabs(x) * DBL_EPSILON / (s * s);
                result->err += result_c.err;
                result->err += DBL_EPSILON * fabs(result->val);
                return; // GSL_SUCCESS;
            }
        }
        else
        {
            result->val = log(y) - 0.5 / x + result_c.val;
            result->err = result_c.err;
            result->err += DBL_EPSILON * fabs(result->val);
            return; // GSL_SUCCESS;
        }
    }
    else
    { // -2 < x < 2
        SFResult result_c;

        if (x < -1.0)
        { // x = -2 + v
            const BS_REAL v  = x + 2.0;
            const BS_REAL t1 = 1.0 / x;
            const BS_REAL t2 = 1.0 / (x + 1.0);
            const BS_REAL t3 = 1.0 / v;
            ChebEvalE(&psi_cs, 2.0 * v - 1.0, &result_c);

            result->val = -(t1 + t2 + t3) + result_c.val;
            result->err = DBL_EPSILON * (fabs(t1) + fabs(x / (t2 * t2)) +
                                         fabs(x / (t3 * t3)));
            result->err += result_c.err;
            result->err += DBL_EPSILON * fabs(result->val);
            return; // GSL_SUCCESS;
        }
        else if (x < 0.0)
        { // x = -1 + v
            const BS_REAL v  = x + 1.0;
            const BS_REAL t1 = 1.0 / x;
            const BS_REAL t2 = 1.0 / v;
            ChebEvalE(&psi_cs, 2.0 * v - 1.0, &result_c);

            result->val = -(t1 + t2) + result_c.val;
            result->err = DBL_EPSILON * (fabs(t1) + fabs(x / (t2 * t2)));
            result->err += result_c.err;
            result->err += DBL_EPSILON * fabs(result->val);
            return; // GSL_SUCCESS;
        }
        else if (x < 1.0)
        { // x = v
            const BS_REAL t1 = 1.0 / x;
            ChebEvalE(&psi_cs, 2.0 * x - 1.0, &result_c);

            result->val = -t1 + result_c.val;
            result->err = DBL_EPSILON * t1;
            result->err += result_c.err;
            result->err += DBL_EPSILON * fabs(result->val);
            return; // GSL_SUCCESS;
        }
        else
        { // x = 1 + v
            const BS_REAL v = x - 1.0;
            ChebEvalE(&psi_cs, 2.0 * v - 1.0, result);
            return;
        }
    }
}
*/

// Evaluation of Psi (Digamma) function (only result)
/*
KOKKOS_INLINE_FUNCTION
BS_REAL SFPsi(const BS_REAL x)
{
    SFResult result;
    SFPsiOutput(x, &result);
    return result.val;
}
*/

/*===========================================================================*/

// fermi_integrals.c

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = -9/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_m92(const BS_REAL x)
{
    const BS_REAL factor = -2. / 7.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[12] = {
        0.270088205852269109,   98457.373700674016,  7895.1716828987114,
        1170.39154599810719,    84.347808600628057,  1.94779771326323869,
        0.00062804476066590963, 32220.8981607867840, 20608.4545396913128,
        5234.8569367171782,     658.75071261067250,  40.9380138615149699};
    constexpr BS_REAL B[21] = {
        -1.8432713155273516,    -7.2420131297087391,    -2.2155513134917048,
        0.04304761419364448,    -0.16668060134499780,   0.10403025949642735,
        0.016505424797869711,   -0.0063406328104926592, -0.0009078668865214013,
        0.00024692958093140298, 129.8291336259148,      -1.9134679779378,
        249.6123960853833,      -3.3454980273582,       190.7566926415911,
        -2.2190711173347,       72.29929273955658,      -0.66382343556531,
        13.551723146833901,     -0.0757657170023758,    1.82715356570827096};
    constexpr BS_REAL C[19] = {
        27.9376225634323064,   17.6954702922794847,  -18.2160790482833106,
        -2.57045608672749918,  2.20866272949862340,  -0.104202929779896185,
        -0.177712759450918142, 0.055743070525419128, 0.0059642251814009821,
        599.97601432112870,    -64.538219224509587,  1048.07497435809354,
        -104.535376308139931,  694.31701140352634,   -61.402988847817226,
        207.271703307329581,   -14.7241003912299059, 23.5639226214135345,
        0.557104412737456300};
    constexpr BS_REAL D[20] = {
        0.105736863491337914,    0.182686468618144972,   -0.0830333880807117100,
        0.0123471104343183151,   0.0213438566082710205,  -0.0431206240183745280,
        0.0311190249582737067,   -0.0115670162945002128, 0.00222720333535487902,
        0.000179769912955733772, 11.4369029896970731,    29.5790675947434036,
        55.2958932613908644,     59.7192644653455683,    51.8646589384828875,
        26.5061101661226671,     11.1515306816966566,    -0.140562900980389991,
        -0.566848363630388431,   3.79898975457237574};
    constexpr BS_REAL E[21] = {
        0.0131615741211877996,   0.0630819955070757629,
        0.0347025335431948214,   0.00659613983600166143,
        0.0239129779624141424,   0.000528952248906432005,
        0.00411353218448862483,  -0.00119448333536781374,
        0.000497581590580634773, -0.0000934288845342437109,
        7.29489294899966855e-6,  12.2828253912160684,
        63.2143680817101004,     155.880993448869246,
        239.067686638304706,     254.656949674257383,
        200.134181826630681,     120.784221710133443,
        56.5037249583191013,     20.2196452284175417,
        5.22412387722958575};
    constexpr BS_REAL F[21] = {
        0.000243271785585373458,  0.000534699823019827761,
        0.000824373965146324914,  0.00106321349769688454,
        0.000488890808530810197,  0.000510404811590977338,
        0.0000543711407291589868, 0.0000597066095966549458,
        -9.52409541937957252e-6,  1.26499979908255436e-6,
        8.23879904655758674e-8,   1.89001039755418574,
        12.0291520111708943,      36.3016148236874403,
        69.0194326307265861,      92.5327210822279397,
        92.3443444983825834,      70.3437819943334403,
        41.0546070181362401,      17.9813587502434782,
        5.52467389311001165};
    constexpr BS_REAL G[19] = {
        5.76414987722574846e-6,   0.0000175838621123388022,
        0.0000306497936833891724, 0.0000245465851262614617,
        9.68986921084327441e-6,   -9.60429091051474115e-6,
        -2.71258276002580026e-6,  2.14306002093744082e-7,
        -2.02878733208210684e-8,  1.11972866705178940e-9,
        0.672666709110930353,     4.51114986457223474,
        14.2130577502093159,      27.0924755033205922,
        33.3161322251216687,      25.3582658290231325,
        8.61719049364763904,      -2.98613950390747585,
        -3.72974082078125757};
    constexpr BS_REAL H[12] = {1,
                               20360.4093941608469,
                               48690.0157659639679,
                               56148.8553220593957,
                               16506.5714937251887,
                               1066.30280063583769,
                               3.73715435019239458,
                               765.319028032818290,
                               1851.12566867527470,
                               2160.21729234460541,
                               677.261952266814024,
                               56.0179124057448609};
    constexpr BS_REAL half  = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex *
             (A[0] -
              ex *
                  (A[1] +
                   t * (A[2] +
                        t * (A[3] + t * (A[4] + t * (A[5] + t * A[6]))))) /
                  (A[7] +
                   t * (A[8] + t * (A[9] + t * (A[10] + t * (A[11] + t))))));
    }
    else if (x < zero)
    {
        s = -half * x;

        fd = (B[0] +
              s * (B[1] +
                   s * (B[2] +
                        s * (B[3] +
                             s * (B[4] +
                                  s * (B[5] +
                                       s * (B[6] +
                                            s * (B[7] +
                                                 s * (B[8] + s * B[9]))))))))) /
             (B[10] +
              s * (B[11] +
                   s * (B[12] +
                        s * (B[13] +
                             s * (B[14] +
                                  s * (B[15] +
                                       s * (B[16] +
                                            s * (B[17] +
                                                 s * (B[18] +
                                                      s * (B[19] + s)))))))))) *
             (x + B[20]); // care for zero point
    }
    else if (x < two)
    {
        t = half * x;

        fd =
            (C[0] +
             t * (C[1] +
                  t * (C[2] +
                       t * (C[3] +
                            t * (C[4] +
                                 t * (C[5] +
                                      t * (C[6] + t * (C[7] - t * C[8])))))))) /
            (C[9] + t * (C[10] +
                         t * (C[11] +
                              t * (C[12] +
                                   t * (C[13] +
                                        t * (C[14] +
                                             t * (C[15] +
                                                  t * (C[16] +
                                                       t * (C[17] - t))))))))) *
            (x - C[18]); // care for zero point
    }
    else if (x < five)
    {
        t = one_third * (x - two);

        fd =
            -(D[0] +
              t * (D[1] +
                   t * (D[2] +
                        t * (D[3] +
                             t * (D[4] +
                                  t * (D[5] +
                                       t * (D[6] +
                                            t * (D[7] +
                                                 t * (D[8] - t * D[9]))))))))) /
            (D[10] +
             t * (D[11] +
                  t * (D[12] +
                       t * (D[13] +
                            t * (D[14] +
                                 t * (D[15] +
                                      t * (D[16] +
                                           t * (D[17] +
                                                t * (D[18] - t))))))))) *
            (x - D[19]); // care for zero point
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;

        fd =
            -(E[0] +
              t * (E[1] +
                   t * (E[2] +
                        t * (E[3] +
                             t * (E[4] +
                                  t * (E[5] +
                                       t * (E[6] +
                                            t * (E[7] +
                                                 t * (E[8] +
                                                      t * (E[9] +
                                                           t * E[10])))))))))) /
            (E[11] +
             t * (E[12] +
                  t * (E[13] +
                       t * (E[14] +
                            t * (E[15] +
                                 t * (E[16] +
                                      t * (E[17] +
                                           t * (E[18] +
                                                t * (E[19] +
                                                     t * (E[20] + t))))))))));
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;

        fd =
            -(F[0] +
              t * (F[1] +
                   t * (F[2] +
                        t * (F[3] +
                             t * (F[4] +
                                  t * (F[5] +
                                       t * (F[6] +
                                            t * (F[7] +
                                                 t * (F[8] +
                                                      t * (F[9] -
                                                           t * F[10])))))))))) /
            (F[11] +
             t * (F[12] +
                  t * (F[13] +
                       t * (F[14] +
                            t * (F[15] +
                                 t * (F[16] +
                                      t * (F[17] +
                                           t * (F[18] +
                                                t * (F[19] +
                                                     t * (F[20] + t))))))))));
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;

        fd =
            -(G[0] +
              t * (G[1] +
                   t * (G[2] +
                        t * (G[3] +
                             t * (G[4] +
                                  t * (G[5] +
                                       t * (G[6] +
                                            t * (G[7] +
                                                 t * (G[8] + t * G[9]))))))))) /
            (G[10] +
             t * (G[11] +
                  t * (G[12] +
                       t * (G[13] +
                            t * (G[14] +
                                 t * (G[15] +
                                      t * (G[16] +
                                           t * (G[17] + t * (G[18] - t)))))))));
    }
    else
    {
        w  = one / (x * x);
        s  = one - onethousandsixhundred * w;
        fd = factor / (sqrt(x) * x * x * x) *
             (H[0] +
              w *
                  (H[1] +
                   s * (H[2] +
                        s * (H[3] + s * (H[4] + s * (H[5] + s * H[6]))))) /
                  (H[7] +
                   s * (H[8] + s * (H[9] + s * (H[10] + s * (H[11] + s))))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = -7/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_m72(const BS_REAL x)
{
    const BS_REAL factor = -2. / 5.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[10] = {-0.945308720482941881, 22521.6817696480830,
                               3006.50943115266211,   223.777970320362539,
                               7.1096231598119406,    0.00186000064211079826,
                               4211.64894546505314,   2132.92301652550982,
                               400.935813755588393,   33.0304914363496116};
    constexpr BS_REAL B[18] = {
        5.49736264373955718,      3.93868182790099018,
        0.00680535603918657089,   -0.0610218211998865281,
        -0.195616572333524455,    0.0172573301744899007,
        0.0148839062617757320,    -0.00259857880173559174,
        0.0000170523346905109458, 48.4238197836838959,
        8.18206322265668686,      73.9370055146725832,
        10.7170918787834126,      42.0907078084017798,
        4.72976153773643018,      10.5965223842451046,
        0.704336166377684061,     0.731435761340666539};
    constexpr BS_REAL C[19] = {
        143.837264704139569,   322.683997310438763,  -209.374575032955500,
        -54.681932464048843,   -30.3617985514093593, -5.7605448223972412,
        3.65745825030136209,   0.145440718467921945, -0.340001125998911237,
        0.0446577756904410773, 1732.20737492411446,  97.952930820143076,
        2660.50150777614871,   115.730278615026183,  1524.22026694958607,
        40.4116469657132500,   385.907052068138736,  1.45068964456455148,
        36.4544527882983039};
    constexpr BS_REAL D[20] = {
        66.1873851488861175,    57.9197865216898287,   -10.6658833177171402,
        10.2066946495410176,    1.19595124198914370,   -1.81286879962085966,
        1.52800227451442796,    -0.546498161708704000, 0.0988077919447287808,
        0.00745279173673929794, 1525.13818882033488,   3966.07893650424194,
        6837.27734004268945,    7086.45653543153192,   5612.71274301484627,
        2928.56588821761611,    1256.17847516062581,   310.558929999645079,
        89.3091775022809646,    2.59355543171509487};
    constexpr BS_REAL E[19] = {
        0.184179377827614791,    0.516542442597467614,
        0.292799259677321999,    0.166303092691678597,
        0.258615034442190300,    -0.0303157249245734253,
        0.0373447714795897992,   -0.00332926057033243781,
        0.000563408857595472728, 0.0000334037183021934960,
        18.2615126816087432,     85.1686926876657081,
        181.378990080239764,     232.528454347497368,
        200.564064511137621,     125.206874911765040,
        58.6883968744415564,     20.6768210567814309,
        5.02620698086580395};
    constexpr BS_REAL F[21] = {0.00521027824483942692,  0.0159223834333424812,
                               0.0290786393625134859,   0.0350554456604475612,
                               0.0276022017688728517,   0.0190056432579866130,
                               0.00725681379911524042,  0.00278845414163336780,
                               0.000285395695659601419, -9.95541819385254380e-6,
                               2.96666654254966665e-7,  3.40047012308625195,
                               20.3896764592442369,     58.0983301256082527,
                               104.466101980605120,     132.574327309248600,
                               125.279491711147983,     90.3637198347191251,
                               49.8942487159252830,     20.6309250519127594,
                               5.96507377885482873};
    constexpr BS_REAL G[19] = {
        0.000222755850099679261, 0.00115329637370955141,
        0.00290802444182764526,  0.00445062601336644735,
        0.00436605816073874863,  0.00260703930450196953,
        0.000708204675930321953, 0.0000387281377930932522,
        -1.07635110532534644e-6, 3.45876932703497152e-8,
        0.958823449511661332,    7.43981964725756497,
        27.2010614567095770,     61.4593095996041459,
        94.2239246717949674,     100.959076669111257,
        74.4305466149867624,     35.3635458846636560,
        9.43423871492993288};
    constexpr BS_REAL H[10] = {1,
                               114389.966467649306,
                               109888.157408212290,
                               21209.9082486576965,
                               896.625771647944984,
                               2.27005555539139462,
                               7803.86522209355981,
                               7642.62657004470466,
                               1584.72834559231955,
                               86.0988387732258338};
    constexpr BS_REAL half  = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex *
             (A[0] +
              ex * (A[1] + t * (A[2] + t * (A[3] + t * (A[4] - t * A[5])))) /
                  (A[6] + t * (A[7] + t * (A[8] + t * (A[9] + t)))));
    }
    else if (x < zero)
    {
        s = -half * x;

        fd =
            (B[0] +
             s * (B[1] +
                  s * (B[2] +
                       s * (B[3] +
                            s * (B[4] +
                                 s * (B[5] +
                                      s * (B[6] + s * (B[7] + s * B[8])))))))) /
            (B[9] +
             s * (B[10] +
                  s * (B[11] +
                       s * (B[12] +
                            s * (B[13] +
                                 s * (B[14] +
                                      s * (B[15] + s * (B[16] + s)))))))) *
            (x + B[17]); // care for zero point
    }
    else if (x < two)
    {
        t = half * x;

        fd =
            (C[0] +
             t * (C[1] +
                  t * (C[2] +
                       t * (C[3] +
                            t * (C[4] +
                                 t * (C[5] +
                                      t * (C[6] +
                                           t * (C[7] +
                                                t * (C[8] + t * C[9]))))))))) /
            (C[10] +
             t * (C[11] +
                  t * (C[12] +
                       t * (C[13] +
                            t * (C[14] +
                                 t * (C[15] +
                                      t * (C[16] +
                                           t * (C[17] + t * (C[18] - t)))))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);

        fd =
            -(D[0] +
              t * (D[1] +
                   t * (D[2] +
                        t * (D[3] +
                             t * (D[4] +
                                  t * (D[5] +
                                       t * (D[6] +
                                            t * (D[7] +
                                                 t * (D[8] - t * D[9]))))))))) /
            (D[10] +
             t * (D[11] +
                  t * (D[12] +
                       t * (D[13] +
                            t * (D[14] +
                                 t * (D[15] +
                                      t * (D[16] +
                                           t * (D[17] +
                                                t * (D[18] - t))))))))) *
            (x - D[19]); // care for zero point
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;

        fd =
            -(E[0] +
              t * (E[1] +
                   t * (E[2] +
                        t * (E[3] +
                             t * (E[4] +
                                  t * (E[5] +
                                       t * (E[6] +
                                            t * (E[7] +
                                                 t * (E[8] - t * E[9]))))))))) /
            (E[10] +
             t * (E[11] +
                  t * (E[12] +
                       t * (E[13] +
                            t * (E[14] +
                                 t * (E[15] +
                                      t * (E[16] +
                                           t * (E[17] + t * (E[18] + t)))))))));
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;

        fd =
            -(F[0] +
              t * (F[1] +
                   t * (F[2] +
                        t * (F[3] +
                             t * (F[4] +
                                  t * (F[5] +
                                       t * (F[6] +
                                            t * (F[7] +
                                                 t * (F[8] +
                                                      t * (F[9] +
                                                           t * F[10])))))))))) /
            (F[11] +
             t * (F[12] +
                  t * (F[13] +
                       t * (F[14] +
                            t * (F[15] +
                                 t * (F[16] +
                                      t * (F[17] +
                                           t * (F[18] +
                                                t * (F[19] +
                                                     t * (F[20] + t))))))))));
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;

        fd =
            -(G[0] +
              t * (G[1] +
                   t * (G[2] +
                        t * (G[3] +
                             t * (G[4] +
                                  t * (G[5] +
                                       t * (G[6] +
                                            t * (G[7] +
                                                 t * (G[8] + t * G[9]))))))))) /
            (G[10] +
             t * (G[11] +
                  t * (G[12] +
                       t * (G[13] +
                            t * (G[14] +
                                 t * (G[15] +
                                      t * (G[16] +
                                           t * (G[17] + t * (G[18] + t)))))))));
    }
    else
    {
        w  = one / (x * x);
        s  = one - onethousandsixhundred * w;
        fd = factor / (sqrt(x) * x * x) *
             (H[0] +
              w * (H[1] + s * (H[2] + s * (H[3] + s * (H[4] + s * H[5])))) /
                  (H[6] + s * (H[7] + s * (H[8] + s * (H[9] + s)))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k=-5/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_m52(const BS_REAL x)
{
    const BS_REAL factor = -2. / 3.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[10] = {2.36327180120735470, 55891.891151428516,
                               10900.4507673238614, 731.61871282829311,
                               17.3576202287841595, 0.00111269690244199535,
                               8361.6144419623308,  3709.66537451151168,
                               598.60597628606205,  41.1772384480908825};
    constexpr BS_REAL B[17] = {
        46.7273884990861909,   65.1460094720913271,   -23.3560330148569366,
        16.3832748811699383,   -6.19159773777140396,  0.799041627235852088,
        -0.136702158227027886, 0.0592594484852497570, 0.00813537888028213710,
        166.600811683081505,   -14.3475443594245037,  207.623210767457661,
        -16.0651771000195284,  91.9265825199149767,   -5.91543082309762067,
        16.6946542439476255,   -0.714551880206870314};
    constexpr BS_REAL C[17] = {
        55.096201117757232,   14.3294333176185961,    8.7057675255700556,
        2.85717540287481278,  0.153440611811510950,   0.111339313377987928,
        0.052394928476702755, 0.00308189450941126788, 217.840649625384406,
        -13.7539547994891024, 264.475360271691825,    -16.3936568402444332,
        112.575929318729487,  -6.5735021347402422,    19.1324946771515682,
        -0.89572129073702419, 1.10894923342229868};
    constexpr BS_REAL D[17] = {
        4.91753613607304828,  20.4599278677959969,   14.7299603843029675,
        6.07034924563989291,  2.84749597355562421,   1.42067667370006629,
        0.194025010869198784, 0.0296942990922311325, 0.0000687311336894144097,
        39.1942895728295075,  102.720679402185990,   163.814228574815100,
        160.042363559709804,  112.736106158227044,   53.5881251269350506,
        19.3013585080264391,  4.60099930165203599};
    constexpr BS_REAL E[17] = {
        2.27928384600199307,  4.20924981693908136,   3.71772880855517432,
        3.45602273475008183,  1.35886950179549082,   0.536408369715153184,
        0.133773439887674841, 0.0125358162552656152, 0.000258659892241001720,
        29.5538409791709387,  102.889070328418668,   171.483234327128671,
        174.069526985568413,  121.638517061339171,   61.2471517251153565,
        22.6513730231817499,  5.63585485331701445};
    constexpr BS_REAL F[19] = {
        0.138865047109769882,     0.495879942014454703,  0.975841030964797014,
        1.19756546203286062,      1.04262591465861665,   0.636815947333856586,
        0.267431185397286764,     0.0680541575349993884, 0.00339987875007404643,
        0.0000596473600392971721, 6.10307238045677639,   32.0683844393622288,
        81.7709446404731603,      131.907213444156949,   148.807381678792024,
        122.126271732385001,      73.4713378422568434,   31.3211204388750577,
        8.63324324695581677};
    constexpr BS_REAL G[17] = {
        0.294313222065093082,  1.31885924825320721,     2.84488318326892294,
        3.47529649317584998,   2.44599317338319296,     0.711594188260176831,
        0.0568475228400688404, 0.000310240841471278615, 7.25941094961112814e-6,
        38.8567648628867585,   233.714122704911173,     657.094399788013655,
        1097.58958431940912,   1157.52830939122072,     742.054932624983009,
        251.716945386093114,   34.4183681886106367};
    constexpr BS_REAL H[9] = {1,
                              563249.531577994933,
                              348717.634501234702,
                              40765.2886083808789,
                              865.059113888543852,
                              90262.6676538587293,
                              56944.6907876652433,
                              7175.67113380361324,
                              207.376531486457753};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex *
             (A[0] -
              ex * (A[1] + t * (A[2] + t * (A[3] + t * (A[4] + t * A[5])))) /
                  (A[6] + t * (A[7] + t * (A[8] + t * (A[9] + t)))));
    }
    else if (x < zero)
    {
        s = -half * x;

        fd =
            (B[0] +
             s * (B[1] +
                  s * (B[2] +
                       s * (B[3] +
                            s * (B[4] +
                                 s * (B[5] +
                                      s * (B[6] + s * (B[7] - s * B[8])))))))) /
            (B[9] +
             s * (B[10] +
                  s * (B[11] +
                       s * (B[12] +
                            s * (B[13] +
                                 s * (B[14] +
                                      s * (B[15] + s * (B[16] + s))))))));
    }
    else if (x < two)
    {
        t = half * x;

        fd = -(C[0] +
               t * (C[1] +
                    t * (C[2] +
                         t * (C[3] +
                              t * (C[4] +
                                   t * (C[5] + t * (C[6] - t * C[7]))))))) /
             (C[8] +
              t * (C[9] +
                   t * (C[10] +
                        t * (C[11] +
                             t * (C[12] +
                                  t * (C[13] +
                                       t * (C[14] + t * (C[15] + t)))))))) *
             (x - C[16]); // care for zero point
    }
    else if (x < five)
    {
        t = one_third * (x - two);

        fd = -(D[0] +
               t * (D[1] +
                    t * (D[2] +
                         t * (D[3] +
                              t * (D[4] +
                                   t * (D[5] +
                                        t * (D[6] +
                                             t * (D[7] + t * D[8])))))))) /
             (D[9] +
              t * (D[10] +
                   t * (D[11] +
                        t * (D[12] +
                             t * (D[13] +
                                  t * (D[14] +
                                       t * (D[15] + t * (D[16] + t))))))));
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;

        fd = -(E[0] +
               t * (E[1] +
                    t * (E[2] +
                         t * (E[3] +
                              t * (E[4] +
                                   t * (E[5] +
                                        t * (E[6] +
                                             t * (E[7] - t * E[8])))))))) /
             (E[9] +
              t * (E[10] +
                   t * (E[11] +
                        t * (E[12] +
                             t * (E[13] +
                                  t * (E[14] +
                                       t * (E[15] + t * (E[16] + t))))))));
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;

        fd =
            -(F[0] +
              t * (F[1] +
                   t * (F[2] +
                        t * (F[3] +
                             t * (F[4] +
                                  t * (F[5] +
                                       t * (F[6] +
                                            t * (F[7] +
                                                 t * (F[8] - t * F[9]))))))))) /
            (F[10] +
             t * (F[11] +
                  t * (F[12] +
                       t * (F[13] +
                            t * (F[14] +
                                 t * (F[15] +
                                      t * (F[16] +
                                           t * (F[17] + t * (F[18] + t)))))))));
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;

        fd = -(G[0] +
               t * (G[1] +
                    t * (G[2] +
                         t * (G[3] +
                              t * (G[4] +
                                   t * (G[5] +
                                        t * (G[6] +
                                             t * (G[7] + t * G[8])))))))) /
             (G[9] +
              t * (G[10] +
                   t * (G[11] +
                        t * (G[12] +
                             t * (G[13] +
                                  t * (G[14] +
                                       t * (G[15] + t * (G[16] + t))))))));
    }
    else
    {
        w = one / (x * x);
        s = one - onethousandsixhundred * w;

        fd = factor / (sqrt(x) * x) *
             (H[0] + w * (H[1] + s * (H[2] + s * (H[3] + s * H[4]))) /
                         (H[5] + s * (H[6] + s * (H[7] + s * (H[8] + s)))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = -3/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_m32(const BS_REAL x)
{
    const BS_REAL factor = -2.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[10] = {-3.54490770181103205, 82737.595643818605,
                               18481.5553495836940,  1272.73919064487495,
                               26.3420403338352574,  0.00110648970639283347,
                               16503.7625405383183,  6422.0552658801394,
                               890.85389683932154,   51.251447078851450};
    constexpr BS_REAL B[17] = {
        946.638483706348559, 76.3328330396778450,  62.7809183134124193,
        83.8442376534073219, 23.2285755924515097,  3.21516808559640925,
        1.58754232369392539, 0.687397326417193593, 0.111510355441975495,
        889.4123665319664,   126.7054690302768,    881.4713137175090,
        108.2557767973694,   289.38131234794585,   27.75902071820822,
        34.252606975067480,  1.9592981990370705};
    constexpr BS_REAL C[17] = {
        754.61690882095729,  565.56180911009650,  494.901267018948095,
        267.922900418996927, 110.418683240337860, 39.4050164908951420,
        10.8654460206463482, 2.11194887477009033, 0.246843599687496060,
        560.03894899770103,  70.007586553114572,  582.42052644718871,
        56.181678606544951,  205.248662395572799, 12.5169006932790528,
        27.2916998671096202, 0.53299717876883183};
    constexpr BS_REAL D[17] = {
        526.022770226139287, 631.116211478274904,  516.367876532501329,
        267.894697896892166, 91.3331816844847913,  17.5723541971644845,
        1.46434478819185576, 1.29615441010250662,  0.223495452221465265,
        354.867400305615304, 560.931137013002977,  666.070260050472570,
        363.745894096653220, 172.272943258816724,  23.7751062504377332,
        12.5916012142616255, -0.888604976123420661};
    constexpr BS_REAL E[14] = {
        18.0110784494455205,  36.1225408181257913,    38.4464752521373310,
        24.1477896166966673,  9.27772356782901602,    2.49074754470533706,
        0.163824586249464178, 0.00329391807590771789, 18.8976860386360201,
        49.3696375710309920,  60.9273314194720251,    43.6334649971575003,
        20.6568810936423065,  6.11094689399482273};
    constexpr BS_REAL F[16] = {
        4.10698092142661427, 17.1412152818912658,   32.6347877674122945,
        36.6653101837618939, 25.9424894559624544,   11.2179995003884922,
        2.30099511642112478, 0.0928307248942099967, 0.00146397877054988411,
        6.40341731836622598, 30.1333068545276116,   64.0494725642004179,
        80.5635003792282196, 64.9297873014508805,   33.3013900893183129,
        9.61549304470339929};
    constexpr BS_REAL G[16] = {
        95.2141371910496454,  420.050572604265456,   797.778374374075796,
        750.378359146985564,  324.818150247463736,   50.3115388695905757,
        0.372431961605507103, -0.103162211894757911, 0.00191752611445211151,
        212.232981736099697,  1043.79079070035083,   2224.50099218470684,
        2464.84669868672670,  1392.55318009810070,   346.597189642259199,
        22.7314613168652593};
    constexpr BS_REAL H[8] = {1,
                              12264.3569103180524,
                              3204.34872454052352,
                              140.119604748253961,
                              0.523918919699235590,
                              9877.87829948067200,
                              2644.71979353906092,
                              128.863768007644572};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex *
             (A[0] +
              ex * (A[1] + t * (A[2] + t * (A[3] + t * (A[4] - t * A[5])))) /
                  (A[6] + t * (A[7] + t * (A[8] + t * (A[9] + t)))));
    }
    else if (x < zero)
    {
        s = -half * x;
        t = one - s;

        fd = -(B[0] +
               t * (B[1] +
                    t * (B[2] +
                         t * (B[3] +
                              t * (B[4] +
                                   t * (B[5] +
                                        t * (B[6] +
                                             t * (B[7] + t * B[8])))))))) /
             (B[9] +
              s * (B[10] +
                   s * (B[11] +
                        s * (B[12] +
                             s * (B[13] +
                                  s * (B[14] +
                                       s * (B[15] + s * (B[16] + s))))))));
    }
    else if (x < two)
    {
        t = half * x;

        fd = -(C[0] +
               t * (C[1] +
                    t * (C[2] +
                         t * (C[3] +
                              t * (C[4] +
                                   t * (C[5] +
                                        t * (C[6] +
                                             t * (C[7] + t * C[8])))))))) /
             (C[9] +
              t * (C[10] +
                   t * (C[11] +
                        t * (C[12] +
                             t * (C[13] +
                                  t * (C[14] +
                                       t * (C[15] + t * (C[16] + t))))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);

        fd = -(D[0] +
               t * (D[1] +
                    t * (D[2] +
                         t * (D[3] +
                              t * (D[4] +
                                   t * (D[5] +
                                        t * (D[6] +
                                             t * (D[7] + t * D[8])))))))) /
             (D[9] +
              t * (D[10] +
                   t * (D[11] +
                        t * (D[12] +
                             t * (D[13] +
                                  t * (D[14] +
                                       t * (D[15] + t * (D[16] + t))))))));
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;

        fd = -(E[0] +
               t * (E[1] +
                    t * (E[2] +
                         t * (E[3] +
                              t * (E[4] +
                                   t * (E[5] + t * (E[6] - t * E[7]))))))) /
             (E[8] +
              t * (E[9] +
                   t * (E[10] + t * (E[11] + t * (E[12] + t * (E[13] + t))))));
    }
    else if (x < twenty)
    {
        t  = one_tenth * x - one;
        fd = -(F[0] +
               t * (F[1] +
                    t * (F[2] +
                         t * (F[3] +
                              t * (F[4] +
                                   t * (F[5] +
                                        t * (F[6] +
                                             t * (F[7] - t * F[8])))))))) /
             (F[9] +
              t * (F[10] +
                   t * (F[11] +
                        t * (F[12] +
                             t * (F[13] + t * (F[14] + t * (F[15] + t)))))));
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;

        fd = -(G[0] +
               t * (G[1] +
                    t * (G[2] +
                         t * (G[3] +
                              t * (G[4] +
                                   t * (G[5] +
                                        t * (G[6] +
                                             t * (G[7] + t * G[8])))))))) /
             (G[9] +
              t * (G[10] +
                   t * (G[11] +
                        t * (G[12] +
                             t * (G[13] + t * (G[14] + t * (G[15] - t)))))));
    }
    else
    {
        w = one / (x * x);
        s = one - onethousandsixhundred * w;

        fd = factor / sqrt(x) *
             (H[0] + w * (H[1] + s * (H[2] + s * (H[3] + s * H[4]))) /
                         (H[5] + s * (H[6] + s * (H[7] + s))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = -1/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_m12(const BS_REAL x)
{
    const BS_REAL factor = 2.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[10] = {1.77245385090551603, 40641.4537510284430,
                               9395.7080940846442,  649.96168315267301,
                               12.7972295804758967, 0.00153864350767585460,
                               32427.1884765292940, 11079.9205661274782,
                               1322.96627001478859, 63.738361029333467};
    constexpr BS_REAL B[16] = {
        272.770092131932696,  30.8845653844682850,   -6.43537632380366113,
        14.8747473098217879,  4.86928862842142635,   -1.53265834550673654,
        -1.02698898315597491, -0.177686820928605932, 0.00377141325509246441,
        293.075378187667857,  305.818162686270816,   299.962395449297620,
        207.640834087494249,  92.0384803181851755,   37.0164914112791209,
        7.88500950271420583};
    constexpr BS_REAL C[15] = {
        3531.50360568243046, 6077.5339658420037,  6199.7700433981326,
        4412.78701919567594, 2252.27343092810898, 811.84098649224085,
        191.836401053637121, 23.2881838959183802, 3293.83702584796268,
        1528.97474029789098, 2568.48562814986046, 925.64264653555825,
        574.23248354035988,  132.803859320667262, 29.8447166552102115};
    constexpr BS_REAL D[15] = {
        4060.70753404118265, 10812.7291333052766, 13897.5649482242583,
        10628.4749852740029, 5107.70670190679021, 1540.84330126003381,
        284.452720112970331, 29.5214417358484151, 1564.58195612633534,
        2825.75172277850406, 3189.16066169981562, 1955.03979069032571,
        828.000333691814748, 181.498111089518376, 32.0352857794803750};
    constexpr BS_REAL E[14] = {
        1198.41719029557508, 3263.51454554908654,  3874.97588471376487,
        2623.13060317199813, 1100.41355637121217,  267.469532490503605,
        25.4207671812718340, 0.389887754234555773, 273.407957792556998,
        595.918318952058643, 605.202452261660849,  343.183302735619981,
        122.187622015695729, 20.9016359079855933};
    constexpr BS_REAL F[15] = {
        9446.00169435237637, 36843.4448474028632, 63710.1115419926191,
        62985.2197361074768, 37634.5231395700921, 12810.9898627807754,
        1981.56896138920963, 81.4930171897667580, 1500.04697810133666,
        5086.91381052794059, 7730.01593747621895, 6640.83376239360596,
        3338.99590300826393, 860.499043886802984, 78.8565824186926692};
    constexpr BS_REAL G[15] = {
        22977.9657855367223, 123416.616813887781, 261153.765172355107,
        274618.894514095795, 149710.718389924860, 40129.3371700184546,
        4470.46495881415076, 132.684346831002976, 2571.68842525335676,
        12521.4982290775358, 23268.1574325055341, 20477.2320119758141,
        8726.52577962268114, 1647.42896896769909, 106.475275142076623};
    constexpr BS_REAL H[7] = {1,
                              0.411233516712009968,
                              0.00110980410034088951,
                              0.0000113689298990173683,
                              2.56931790679436797e-7,
                              9.97897786755446178e-9,
                              8.67667698791108582e-10};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex *
             (A[0] -
              ex * (A[1] + t * (A[2] + t * (A[3] + t * (A[4] + t * A[5])))) /
                  (A[6] + t * (A[7] + t * (A[8] + t * (A[9] + t)))));
    }
    else if (x < zero)
    {
        s = -half * x;
        t = one - s;

        fd =
            (B[0] +
             t * (B[1] +
                  t * (B[2] +
                       t * (B[3] +
                            t * (B[4] +
                                 t * (B[5] +
                                      t * (B[6] + t * (B[7] - t * B[8])))))))) /
            (B[9] +
             s * (B[10] +
                  s * (B[11] +
                       s * (B[12] +
                            s * (B[13] + s * (B[14] + s * (B[15] + s)))))));
    }
    else if (x < two)
    {
        t = half * x;

        fd = (C[0] +
              t * (C[1] +
                   t * (C[2] +
                        t * (C[3] +
                             t * (C[4] +
                                  t * (C[5] + t * (C[6] + t * C[7]))))))) /
             (C[8] +
              t * (C[9] +
                   t * (C[10] +
                        t * (C[11] +
                             t * (C[12] + t * (C[13] + t * (C[14] + t)))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);

        fd = (D[0] +
              t * (D[1] +
                   t * (D[2] +
                        t * (D[3] +
                             t * (D[4] +
                                  t * (D[5] + t * (D[6] + t * D[7]))))))) /
             (D[8] +
              t * (D[9] +
                   t * (D[10] +
                        t * (D[11] +
                             t * (D[12] + t * (D[13] + t * (D[14] + t)))))));
    }
    else if (x < ten)
    {
        t  = one_fifth * x - one;
        fd = (E[0] +
              t * (E[1] +
                   t * (E[2] +
                        t * (E[3] +
                             t * (E[4] +
                                  t * (E[5] + t * (E[6] + t * E[7]))))))) /
             (E[8] +
              t * (E[9] +
                   t * (E[10] + t * (E[11] + t * (E[12] + t * (E[13] + t))))));
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;

        fd = (F[0] +
              t * (F[1] +
                   t * (F[2] +
                        t * (F[3] +
                             t * (F[4] +
                                  t * (F[5] + t * (F[6] + t * F[7]))))))) /
             (F[8] +
              t * (F[9] +
                   t * (F[10] +
                        t * (F[11] +
                             t * (F[12] + t * (F[13] + t * (F[14] + t)))))));
    }
    else if (x < forty)
    {
        t  = one_twentyth * x - one;
        fd = (G[0] +
              t * (G[1] +
                   t * (G[2] +
                        t * (G[3] +
                             t * (G[4] +
                                  t * (G[5] + t * (G[6] + t * G[7]))))))) /
             (G[8] +
              t * (G[9] +
                   t * (G[10] +
                        t * (G[11] +
                             t * (G[12] + t * (G[13] + t * (G[14] + t)))))));
    }
    else
    {
        w = one / (x * x);
        t = onethousandsixhundred * w;

        fd = sqrt(x) * factor *
             (H[0] -
              w * (H[1] +
                   t * (H[2] +
                        t * (H[3] + t * (H[4] + t * (H[5] + t * H[6]))))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 0 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_0(const BS_REAL y)
{
    BS_REAL x, ex, t, s;
    BS_REAL fd = 0.;

    x = -fabs(y);

    constexpr BS_REAL A[10] = {1,
                               22696.2126132366633,
                               5222.0667923565138,
                               357.623326425354522,
                               6.9167792879948140,
                               0.00200096064827815813,
                               45392.4252264733267,
                               14539.5980679273792,
                               1611.36476693109675,
                               71.072178562726798};
    constexpr BS_REAL B[16] = {
        159.601717762460980,   23.7193942338278703,    0.377783268730614356,
        10.5181677709577503,   3.78181326142271599,    -0.441998676614933572,
        -0.450072959113928254, -0.0734798777056723512, 0.000915454570009894267,
        284.26032127745967,    315.2592651624449,      310.2713981221035,
        206.21640678892182,    96.77898293084927,      35.456591489081173,
        8.1762315442738975};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex *
             (A[0] -
              ex * (A[1] + t * (A[2] + t * (A[3] + t * (A[4] + t * A[5])))) /
                  (A[6] + t * (A[7] + t * (A[8] + t * (A[9] + t)))));
    }
    else if (x <= zero)
    {
        s = -half * x;
        t = one - s;

        fd =
            (B[0] +
             t * (B[1] +
                  t * (B[2] +
                       t * (B[3] +
                            t * (B[4] +
                                 t * (B[5] +
                                      t * (B[6] + t * (B[7] + t * B[8])))))))) /
            (B[9] +
             s * (B[10] +
                  s * (B[11] +
                       s * (B[12] +
                            s * (B[13] + s * (B[14] + s * (B[15] + s)))))));
    }
    if (y > zero)
    {
        fd = fd + y;
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 1/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p12(const BS_REAL x)
{
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[10] = {0.886226925452758014, 19894.4553386951666,
                               4509.64329955948557,  303.461789035142376,
                               5.7574879114754736,   0.00275088986849762610,
                               63493.915041308052,   19070.1178243603945,
                               1962.19362141235102,  79.250704958640158};
    constexpr BS_REAL B[16] = {
        149.462587768865243,   22.8125889885050154,    -0.629256395534285422,
        9.08120441515995244,   3.35357478401835299,    -0.473677696915555805,
        -0.467190913556185953, -0.0880610317272330793, 0.00262208080491572673,
        269.94660938022644,    343.6419926336247,      323.9049470901941,
        218.89170769294024,    102.31331350098315,     36.319337289702664,
        8.3317401231389461};
    constexpr BS_REAL C[17] = {
        71652.717119215557,  134954.734070223743, 153693.833350315645,
        123247.280745703400, 72886.293647930726,  32081.2499422362952,
        10210.9967337762918, 2152.71110381320778, 232.906588165205042,
        105667.839854298798, 31946.0752989314444, 71158.788776422211,
        15650.8990138187414, 13521.8033657783433, 1646.98258283527892,
        618.90691969249409,  -3.36319591755394735};
    constexpr BS_REAL D[15] = {
        23744.8706993314289,  68257.8589855623002,  89327.4467683334597,
        62766.3415600442563,  20093.6622609901994,  -2213.89084119777949,
        -3901.66057267577389, 948.642895944858861,  9488.61972919565851,
        12514.8125526953073,  9903.44088207450946,  2138.15420910334305,
        -528.394863730838233, -661.033633995449691, -51.4481470250962337};
    constexpr BS_REAL E[16] = {
        311337.452661582536,   1.11267074416648198e6, 1.75638628895671735e6,
        1.59630855803772449e6, 910818.935456183774,   326492.733550701245,
        65507.2624972852908,   4809.45649527286889,   39721.6641625089685,
        86424.7529107662431,   88163.7255252151780,   50615.7363511157353,
        17334.9774805008209,   2712.13170809042550,   82.2205828354629102,
        0.999999999999999877};
    constexpr BS_REAL F[14] = {
        7.26870063003059784e6, 2.79049734854776025e7, 4.42791767759742390e7,
        3.63735017512363365e7, 1.55766342463679795e7, 2.97469357085299505e6,
        154516.447031598403,   340542.544360209743,   805021.468647620047,
        759088.235455002605,   304686.671371640343,   39289.4061400542309,
        582.426138126398363,   11.2728194581586028};
    constexpr BS_REAL G[13] = {
        4.81449797541963104e6, 1.85162850713127602e7, 2.77630967522574435e7,
        2.03275937688070624e7, 7.41578871589369361e6, 1.21193113596189034e6,
        63211.9545144644852,   80492.7765975237449,   189328.678152654840,
        151155.890651482570,   48146.3242253837259,   5407.08878394180588,
        112.195044410775577};
    constexpr BS_REAL H[7] = {0.666666666666666667, 1,
                              8109.79390744477921,  342.069867454704106,
                              1.07141702293504595,  6569.98472532829094,
                              280.706465851683809};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex *
             (A[0] -
              ex * (A[1] + t * (A[2] + t * (A[3] + t * (A[4] + t * A[5])))) /
                  (A[6] + t * (A[7] + t * (A[8] + t * (A[9] + t)))));
    }
    else if (x < zero)
    {
        s = -half * x;
        t = one - s;

        fd =
            (B[0] +
             t * (B[1] +
                  t * (B[2] +
                       t * (B[3] +
                            t * (B[4] +
                                 t * (B[5] +
                                      t * (B[6] + t * (B[7] - t * B[8])))))))) /
            (B[9] +
             s * (B[10] +
                  s * (B[11] +
                       s * (B[12] +
                            s * (B[13] + s * (B[14] + s * (B[15] + s)))))));
    }
    else if (x < two)
    {
        t = half * x;

        fd =
            (C[0] +
             t * (C[1] +
                  t * (C[2] +
                       t * (C[3] +
                            t * (C[4] +
                                 t * (C[5] +
                                      t * (C[6] + t * (C[7] + t * C[8])))))))) /
            (C[9] +
             t * (C[10] +
                  t * (C[11] +
                       t * (C[12] +
                            t * (C[13] +
                                 t * (C[14] +
                                      t * (C[15] + t * (C[16] + t))))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);

        fd = (D[0] +
              t * (D[1] +
                   t * (D[2] +
                        t * (D[3] +
                             t * (D[4] +
                                  t * (D[5] + t * (D[6] - t * D[7]))))))) /
             (D[8] +
              t * (D[9] +
                   t * (D[10] +
                        t * (D[11] +
                             t * (D[12] + t * (D[13] + t * (D[14] + t)))))));
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;

        fd = (E[0] +
              t * (E[1] +
                   t * (E[2] +
                        t * (E[3] +
                             t * (E[4] +
                                  t * (E[5] + t * (E[6] + t * E[7]))))))) /
             (E[8] +
              t * (E[9] +
                   t * (E[10] +
                        t * (E[11] +
                             t * (E[12] + t * (E[13] + t * (E[14] - t))))))) *
             E[15];
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;

        fd = (F[0] +
              t * (F[1] +
                   t * (F[2] +
                        t * (F[3] + t * (F[4] + t * (F[5] + t * F[6])))))) /
             (F[7] +
              t * (F[8] +
                   t * (F[9] +
                        t * (F[10] +
                             t * (F[11] + t * (F[12] + t * (F[13] - t)))))));
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;

        fd = (G[0] +
              t * (G[1] +
                   t * (G[2] +
                        t * (G[3] + t * (G[4] + t * (G[5] + t * G[6])))))) /
             (G[7] +
              t * (G[8] +
                   t * (G[9] + t * (G[10] + t * (G[11] + t * (G[12] - t))))));
    }
    else
    {
        w = one / (x * x);
        s = one - onethousandsixhundred * w;

        fd = x * sqrt(x) * H[0] *
             (H[1] +
              w * (H[2] + s * (H[3] + s * H[4])) / (H[5] + s * (H[6] + s)));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 1 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p1(const BS_REAL y)
{
    BS_REAL x, ex, t, s;
    BS_REAL fd = 0.;

    x = -fabs(y);

    constexpr BS_REAL A[10] = {1,
                               22189.1070807945062,
                               4915.92700908746777,
                               322.901386168881348,
                               5.9897442965804548,
                               0.00397641173774375092,
                               88756.428323178025,
                               25002.3197546553836,
                               2389.06277237306633,
                               88.376214553692756};
    constexpr BS_REAL B[17] = {
        145.488167182330098,  251.392824471576922,   56.6537141912783024,
        17.9918985363509694,  20.1369115558099802,   7.09659390228556164,
        0.199701180197912643, -0.403173132925886253, 0.0792966701498222697,
        606.0757707716040,    374.1806357435014,     252.1367051536344,
        27.2746245830016,     -61.57766112137513,    -53.72117554363975,
        -25.678454878692950,  -7.1995819520154718};
    constexpr BS_REAL C[2] = {1.64493406684822644, 0.5};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex *
             (A[0] -
              ex * (A[1] + t * (A[2] + t * (A[3] + t * (A[4] + t * A[5])))) /
                  (A[6] + t * (A[7] + t * (A[8] + t * (A[9] + t)))));
    }
    else if (x <= zero)
    {
        s = -half * x;
        t = one - s;

        fd =
            (B[0] +
             t * (B[1] +
                  t * (B[2] +
                       t * (B[3] +
                            t * (B[4] +
                                 t * (B[5] +
                                      t * (B[6] + t * (B[7] - t * B[8])))))))) /
            (B[9] +
             s * (B[10] +
                  s * (B[11] +
                       s * (B[12] +
                            s * (B[13] +
                                 s * (B[14] +
                                      s * (B[15] + s * (B[16] - s))))))));
    }
    if (y > zero)
    {
        fd = -fd + C[0] + C[1] * y * y;
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 3/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p32(const BS_REAL x)
{
    const BS_REAL factor = 2. / 5.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[8]  = {1.32934038817913702,   1346.14119046566636,
                               199.946876779712712,   6.5210149677288048,
                               0.0108588591982722183, 5728.3481201778541,
                               1132.17837281710987,   64.805243148002602};
    constexpr BS_REAL B[17] = {
        631.667081787115831,  504.131655805666135,   113.449065431934917,
        56.0939647772947784,  43.3374223200846752,   12.8047010597109577,
        0.219164386586949410, -0.678659552658390139, 0.126533769309899232,
        1180.5112183558028,   1101.0159189871135,    864.4448234404281,
        392.2018227840790,    89.58093202779063,     -9.95066218572899,
        -17.312068771997626,  -6.4116162917822773};
    constexpr BS_REAL C[17] = {
        90122.488639370400,   157095.208147064037,  166879.962599668589,
        125708.597728045460,  69968.278213181390,   29035.3292989055404,
        8736.4439472398517,   1747.16784760309227,  180.132410666734053,
        78176.777123671727,   -1681.44633240543085, 38665.7913035496031,
        -2527.29685826087874, 5062.6683078100048,   -553.21165462054589,
        165.395637981775430,  -18.0295465153725544};
    constexpr BS_REAL D[15] = {
        912944.432058014054,   3.28217091334054338e6, 5.59250227196369585e6,
        5.76136129685687470e6, 3.84331519034749983e6, 1.65284168824947710e6,
        423452.676670436605,   49835.4127241373113,   164873.145721762182,
        257442.511191094986,   225604.160532840884,   99932.1955662320024,
        24761.0878784286761,   1398.26392212830777,   -36.4450237523474167};
    constexpr BS_REAL E[17] = {
        1.88412548327216052e6, 8.08838896259910792e6, 1.56869793001790529e7,
        1.79109792599373447e7, 1.31345142328147214e7, 6.29500412046744325e6,
        1.89326213154091054e6, 312372.643127575407,   18814.7420442630170,
        67768.3347951202583,   147635.914444221358,   151908.303165069423,
        86671.1222110642970,   27855.9481608626219,   3833.22697473114940,
        98.3384567064269554,   0.999999999999999876};
    constexpr BS_REAL F[16] = {
        1.59656593348660977e9,  7.32769737561517060e9, 1.42662658588280191e10,
        1.51238422045169918e10, 9.27233604548095476e9, 3.18834513406577423e9,
        5.36061988605886123e8,  3.03619219668246382e7, 1.18906980815759995e7,
        2.62209219322122975e7,  2.28143701746618729e7, 8.57156701742181874e6,
        1.13860063870524239e6,  27091.7884208687379,   -275.664733379090447,
        0.999999999999999829};
    constexpr BS_REAL G[14] = {
        2.60437581212904589e8, 1.08771546307370080e9, 1.81531350939088943e9,
        1.52833764636304939e9, 6.70684451492750149e8, 1.40870639531414149e8,
        1.04957900377463854e7, 358448.871166784200,   611808.419702466190,
        326307.561591723775,   58407.9904827573816,   2049.50040323021794,
        -39.8767861209088081,  0.999999999999999828};
    constexpr BS_REAL H[6] = {1,
                              6.16739021212286242,
                              0.00111530123694574981,
                              -2.79156524536560815e-6,
                              2.95571462110856359e-8,
                              6.70917556862133933e-10};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex * (A[0] - ex * (A[1] + t * (A[2] + t * (A[3] + t * A[4]))) /
                              (A[5] + t * (A[6] + t * (A[7] + t))));
    }
    else if (x < zero)
    {
        s = -half * x;
        t = one - s;

        fd =
            (B[0] +
             t * (B[1] +
                  t * (B[2] +
                       t * (B[3] +
                            t * (B[4] +
                                 t * (B[5] +
                                      t * (B[6] + t * (B[7] - t * B[8])))))))) /
            (B[9] +
             s * (B[10] +
                  s * (B[11] +
                       s * (B[12] +
                            s * (B[13] +
                                 s * (B[14] +
                                      s * (B[15] + s * (B[16] - s))))))));
    }
    else if (x < two)
    {
        t = half * x;

        fd =
            (C[0] +
             t * (C[1] +
                  t * (C[2] +
                       t * (C[3] +
                            t * (C[4] +
                                 t * (C[5] +
                                      t * (C[6] + t * (C[7] + t * C[8])))))))) /
            (C[9] +
             t * (C[10] +
                  t * (C[11] +
                       t * (C[12] +
                            t * (C[13] +
                                 t * (C[14] +
                                      t * (C[15] + t * (C[16] + t))))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);

        fd = (D[0] +
              t * (D[1] +
                   t * (D[2] +
                        t * (D[3] +
                             t * (D[4] +
                                  t * (D[5] + t * (D[6] + t * D[7]))))))) /
             (D[8] +
              t * (D[9] +
                   t * (D[10] +
                        t * (D[11] +
                             t * (D[12] + t * (D[13] + t * (D[14] + t)))))));
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;

        fd =
            (E[0] +
             t * (E[1] +
                  t * (E[2] +
                       t * (E[3] +
                            t * (E[4] +
                                 t * (E[5] +
                                      t * (E[6] + t * (E[7] + t * E[8])))))))) /
            (E[9] +
             t * (E[10] +
                  t * (E[11] +
                       t * (E[12] +
                            t * (E[13] + t * (E[14] + t * (E[15] - t))))))) *
            E[16]; // correction to remove bias
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;

        fd = (F[0] +
              t * (F[1] +
                   t * (F[2] +
                        t * (F[3] +
                             t * (F[4] +
                                  t * (F[5] + t * (F[6] + t * F[7]))))))) /
             (F[8] +
              t * (F[9] +
                   t * (F[10] +
                        t * (F[11] +
                             t * (F[12] + t * (F[13] + t * (F[14] + t))))))) *
             F[15]; // correction to remove bias
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;

        fd = (G[0] +
              t * (G[1] +
                   t * (G[2] +
                        t * (G[3] + t * (G[4] + t * (G[5] + t * G[6])))))) /
             (G[7] +
              t * (G[8] +
                   t * (G[9] + t * (G[10] + t * (G[11] + t * (G[12] + t)))))) *
             G[13];
    }
    else
    {
        w  = one / (x * x);
        s  = one - onethousandsixhundred * w;
        fd = x * x * sqrt(x) * factor *
             (H[0] +
              w * (H[1] + s * (H[2] + s * (H[3] + s * (H[4] - s * H[5])))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p2(const BS_REAL y)
{
    BS_REAL x, ex, t, s;
    BS_REAL fd = 0.;

    x = -fabs(y);

    constexpr BS_REAL A[8]  = {2,
                               1914.06748184935743,
                               273.085756700981399,
                               8.5861610217850095,
                               0.0161890243763741414,
                               7656.2699273974454,
                               1399.35442210906621,
                               72.929152915475392};
    constexpr BS_REAL B[17] = {
        2711.49678259128843,  1299.85460914884154,  222.606134197895041,
        172.881855215582924,  112.951038040682055,  24.0376482128898634,
        -2.68393549333878715, -2.14077421411719935, 0.326188299771397236,
        2517.1726659917047,   3038.7689794575778,   2541.7823512372631,
        1428.0589853413436,   531.62378035996132,   122.54595216479181,
        8.395768655115050,    -3.9142702096919080};
    constexpr BS_REAL C[2] = {3.28986813369645287, 0.333333333333333333};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex * (A[0] - ex * (A[1] + t * (A[2] + t * (A[3] + t * A[4]))) /
                              (A[5] + t * (A[6] + t * (A[7] + t))));
    }
    else if (x <= zero)
    {
        s = -half * x;
        t = one - s;

        fd =
            (B[0] +
             t * (B[1] +
                  t * (B[2] +
                       t * (B[3] +
                            t * (B[4] +
                                 t * (B[5] +
                                      t * (B[6] + t * (B[7] - t * B[8])))))))) /
            (B[9] +
             s * (B[10] +
                  s * (B[11] +
                       s * (B[12] +
                            s * (B[13] +
                                 s * (B[14] +
                                      s * (B[15] + s * (B[16] - s))))))));
    }
    if (y > zero)
    {
        fd = fd + y * (C[0] + C[1] * y * y);
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 5/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p52(const BS_REAL x)
{
    const BS_REAL factor = 2. / 7.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[8]  = {3.32335097044784255,   3004.41138112148735,
                               409.582975848055860,   12.3422465543050559,
                               0.0252600128832787650, 10227.9400759349389,
                               1729.21912667445507,   82.075224736463695};
    constexpr BS_REAL B[15] = {
        3273.17630052701392, 1698.67676090446982,  498.075924016264615,
        286.399893838802205, 153.449411299292925,  46.3172989461474687,
        7.47269262560340991, 0.512665223813025153, 1934.7654167995427,
        2174.8582387970533,  1803.7716559637348,   970.4325125469826,
        368.29988963760778,  94.85059993048974,    14.845914179579222};
    constexpr BS_REAL C[15] = {
        46746.0667679244140, 67556.152887188131,  58097.414724032516,
        35613.1718933770642, 15859.9235610267359, 5040.9166211327297,
        1064.44778209849372, 117.741363279815688, 13126.266942469915,
        -1313.6945396119670, 2636.3767383046264,  453.2433336267639,
        201.68980019364685,  50.35944043583614,   9.319958361334161};
    constexpr BS_REAL D[15] = {
        2.96446396463084245e6, 1.13985063437515144e7, 2.03551842990364919e7,
        2.18469168036723457e7, 1.52194204967854489e7, 6.89564335059367232e6,
        1.88985082880325956e6, 243639.263893338434,   169113.644493386066,
        249598.556346488591,   162470.264413303189,   48830.6689634566512,
        4080.41491690869679,   -227.369223251313984,  19.1547319985914212};
    constexpr BS_REAL E[16] = {
        5.98759698946959779e7, 2.87159410448355925e8, 5.99247230055005978e8,
        7.05261579500244084e8, 5.04521860843522470e8, 2.18416589668256432e8,
        5.25509373690209805e7, 5.33490422476777307e6, 469653.962099017877,
        972160.483958338969,   696829.387839904810,   200021.997093042021,
        12604.5784254924892,   -464.256227099766735,  26.4866799762414562,
        0.999999999999999786};
    constexpr BS_REAL F[16] = {
        1.11471644837527339e9,  5.74225610798369390e9,  1.26551657157853213e10,
        1.53960861052359422e10, 1.10946777927365321e10, 4.68865671930554474e9,
        1.05763202432846973e9,  9.49124183370767767e7,  1.07734938834750844e6,
        2.05459707311873616e6,  1.39824422108531691e6,  347716.997197363113,
        17894.5194245484999,    -553.162195184268593,   28.1090136251865326,
        0.999999999999999759};
    constexpr BS_REAL G[14] = {
        8.02670881104191218e9,  3.67306963017003546e10, 6.83085243661572356e10,
        6.55107149365528043e10, 3.37406319128317264e10, 8.67814222875818408e9,
        8.43844503352450216e8,  757905.984443885971,    868424.806294231233,
        260670.917865642513,    14550.0472712579662,    -476.120164041067762,
        25.8288614974100332,    0.999999999999999820};
    constexpr BS_REAL H[6] = {1,
                              14.3931730849220041,
                              0.00776862867834285253,
                              3.78966458769690333e-6,
                              2.09248095155530095e-8,
                              3.52097438532254351e-10};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex * (A[0] - ex * (A[1] + t * (A[2] + t * (A[3] + t * A[4]))) /
                              (A[5] + t * (A[6] + t * (A[7] + t))));
    }
    else if (x < zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] +
                             t * (B[4] +
                                  t * (B[5] + t * (B[6] + t * B[7]))))))) /
             (B[8] +
              s * (B[9] +
                   s * (B[10] +
                        s * (B[11] +
                             s * (B[12] + s * (B[13] + s * (B[14] + s)))))));
    }
    else if (x < two)
    {
        t = half * x;
        s = one - t;

        fd = (C[0] +
              t * (C[1] +
                   t * (C[2] +
                        t * (C[3] +
                             t * (C[4] +
                                  t * (C[5] + t * (C[6] + t * C[7]))))))) /
             (C[8] +
              s * (C[9] +
                   s * (C[10] +
                        s * (C[11] +
                             s * (C[12] + s * (C[13] + s * (C[14] + s)))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);

        fd = (D[0] +
              t * (D[1] +
                   t * (D[2] +
                        t * (D[3] +
                             t * (D[4] +
                                  t * (D[5] + t * (D[6] + t * D[7]))))))) /
             (D[8] +
              t * (D[9] +
                   t * (D[10] +
                        t * (D[11] +
                             t * (D[12] + t * (D[13] + t * (D[14] - t)))))));
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;

        fd = (E[0] +
              t * (E[1] +
                   t * (E[2] +
                        t * (E[3] +
                             t * (E[4] +
                                  t * (E[5] + t * (E[6] + t * E[7]))))))) /
             (E[8] +
              t * (E[9] +
                   t * (E[10] +
                        t * (E[11] +
                             t * (E[12] + t * (E[13] + t * (E[14] - t))))))) *
             E[15]; // correction to remove bias
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;

        fd = (F[0] +
              t * (F[1] +
                   t * (F[2] +
                        t * (F[3] +
                             t * (F[4] +
                                  t * (F[5] + t * (F[6] + t * F[7]))))))) /
             (F[8] +
              t * (F[9] +
                   t * (F[10] +
                        t * (F[11] +
                             t * (F[12] + t * (F[13] + t * (F[14] - t))))))) *
             F[15]; // correction to remove bias
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;

        fd = (G[0] +
              t * (G[1] +
                   t * (G[2] +
                        t * (G[3] + t * (G[4] + t * (G[5] + t * G[6])))))) /
             (G[7] +
              t * (G[8] +
                   t * (G[9] + t * (G[10] + t * (G[11] + t * (G[12] - t)))))) *
             G[13]; // correction to remove bias
    }
    else
    {
        w = one / (x * x);
        t = onethousandsixhundred * w;

        fd = x * x * x * sqrt(x) * factor *
             (H[0] +
              w * (H[1] + t * (H[2] + t * (H[3] + t * (H[4] + t * H[5])))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 3 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p3(const BS_REAL y)
{
    BS_REAL x, ex, t, s, y2;
    BS_REAL fd = 0.;

    x = -fabs(y);

    constexpr BS_REAL A[8]  = {6,
                               5121.6401850302408,
                               664.28706260743472,
                               19.0856927562699544,
                               0.0410982603688952131,
                               13657.7071600806539,
                               2136.54222460571183,
                               92.376788603062645};
    constexpr BS_REAL B[15] = {
        7881.24597452900838, 4323.07526636309661, 1260.13125873282465,
        653.359212389160499, 354.630774329461644, 113.373708671587772,
        19.9559488532742796, 1.59407954898394322, 2570.7250703533430,
        2972.7443644211129,  2393.9995533270879,  1259.0724833462608,
        459.86413596901097,  112.60906419590854,  16.468882811659000};
    constexpr BS_REAL C[3] = {11.3643939539669510, 4.93480220054467931, 0.25};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex * (A[0] - ex * (A[1] + t * (A[2] + t * (A[3] + t * A[4]))) /
                              (A[5] + t * (A[6] + t * (A[7] + t))));
    }
    else if (x <= zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] +
                             t * (B[4] +
                                  t * (B[5] + t * (B[6] + t * B[7]))))))) /
             (B[8] +
              s * (B[9] +
                   s * (B[10] +
                        s * (B[11] +
                             s * (B[12] + s * (B[13] + s * (B[14] + s)))))));
    }
    if (y > zero)
    {
        y2 = y * y;
        fd = -fd + C[0] + y2 * (C[1] + y2 * C[2]);
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 7/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p72(const BS_REAL x)
{
    const BS_REAL factor = 2. / 9.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[8]  = {11.6317283965674489,  9371.8868848378853,
                               1152.33581795824468,  31.4225568934398918,
                               0.069508514261022902, 18231.3053891121107,
                               2639.60073595942887,  103.985056337236794};
    constexpr BS_REAL B[16] = {
        12412.9547163940264, 6557.39898920194282, 1745.01251050433134,
        883.903030170203178, 485.680209936338180, 150.395267497504100,
        23.6841796284116177, 1.07158463967401877, 0.157435982722068683,
        1990.3886647182248,  2450.5836860178091,  1959.3508463399357,
        1039.6469110978007,  381.70374149314543,  94.69566180593596,
        14.401931782580504};
    constexpr BS_REAL C[15] = {
        138724.961458493789, 169546.564556791848, 125336.240767414531,
        68469.193938175092,  27692.3811560029092, 8054.8076265209314,
        1569.98144698357220, 162.669853437143155, 7624.675802050973,
        2096.5891309081106,  1735.7563869500611,  662.0787128606726,
        218.91333929478294,  55.27329667089387,   9.904579892966869};
    constexpr BS_REAL D[15] = {
        2.20638302461716245e7, 7.19733683266519997e7, 1.11705692645922359e8,
        1.06328515554623708e8, 6.68969329369869135e7, 2.78374130823093789e7,
        7.11915156250133512e6, 870011.051674696596,   311792.108500755594,
        206109.613235058227,   81814.2443130826819,   3513.24062956494342,
        368.155702298381570,   -106.091905490453452,  15.0961925962965432};
    constexpr BS_REAL E[16] = {
        3.0709767203984838e9,  1.3145330402583295e10, 2.5004908399964440e10,
        2.7299220224056906e10, 1.8396311263825673e10, 7.6098692716640943e9,
        1.7746200575457532e9,  1.7740961200347283e8,  4.46811391083410494e6,
        4.62189039302323844e6, 1.60932873023944955e6, 113458.965071032695,
        -4751.36769214919833,  351.424157693641750,   -24.0392624870048634,
        0.999999999999999779};
    constexpr BS_REAL F[15] = {
        4.91615192873221833e9,  2.30941143669878184e10, 4.52380364776557202e10,
        4.68413409944225048e10, 2.65134543816397678e10, 7.35961621586797944e9,
        5.52546586238663078e8,  7.67154234216335832e7,  550749.688306094769,
        352810.964448809157,    12162.0407442978353,    -4414.91208383702609,
        288.821012224214490,    -21.7191096434783647,   0.999999999999999647};
    constexpr BS_REAL G[16] = {
        1.23356074447405045e12, 6.8241092859415704e12, 1.5939237639946117e13,
        2.0277830436278028e13,  1.5074920317497838e13, 6.4878414656043323e12,
        1.4761927839220043e12,  1.3410778049868313e11, 7.2811112110966361e6,
        8.4187326609915276e6,   2.59806064027658564e6, 153974.958146639672,
        -5666.5977917092119,    387.205694265311079,   -25.1500172248070365,
        0.999999999999999750};
    constexpr BS_REAL H[6] = {1,
                              25.9077115528595532,
                              0.0699176581158797143,
                              -0.0000113689722343055157,
                              -2.69205155161558844e-8,
                              2.81838487282327867e-10};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex * (A[0] - ex * (A[1] + t * (A[2] + t * (A[3] + t * A[4]))) /
                              (A[5] + t * (A[6] + t * (A[7] + t))));
    }
    else if (x < zero)
    {
        s = -half * x;
        t = one - s;

        fd =
            (B[0] +
             t * (B[1] +
                  t * (B[2] +
                       t * (B[3] +
                            t * (B[4] +
                                 t * (B[5] +
                                      t * (B[6] + t * (B[7] - t * B[8])))))))) /
            (B[9] +
             s * (B[10] +
                  s * (B[11] +
                       s * (B[12] +
                            s * (B[13] + s * (B[14] + s * (B[15] + s)))))));
    }
    else if (x < two)
    {
        t = half * x;
        s = one - t;

        fd = (C[0] +
              t * (C[1] +
                   t * (C[2] +
                        t * (C[3] +
                             t * (C[4] +
                                  t * (C[5] + t * (C[6] + t * C[7]))))))) /
             (C[8] +
              s * (C[9] +
                   s * (C[10] +
                        s * (C[11] +
                             s * (C[12] + s * (C[13] + s * (C[14] + s)))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);

        fd = (D[0] +
              t * (D[1] +
                   t * (D[2] +
                        t * (D[3] +
                             t * (D[4] +
                                  t * (D[5] + t * (D[6] + t * D[7]))))))) /
             (D[8] +
              t * (D[9] +
                   t * (D[10] +
                        t * (D[11] +
                             t * (D[12] + t * (D[13] + t * (D[14] - t)))))));
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;

        fd = (E[0] +
              t * (E[1] +
                   t * (E[2] +
                        t * (E[3] +
                             t * (E[4] +
                                  t * (E[5] + t * (E[6] + t * E[7]))))))) /
             (E[8] +
              t * (E[9] +
                   t * (E[10] +
                        t * (E[11] +
                             t * (E[12] + t * (E[13] + t * (E[14] + t))))))) *
             E[15]; // correction to remove bias
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;

        fd = (F[0] +
              t * (F[1] +
                   t * (F[2] +
                        t * (F[3] +
                             t * (F[4] +
                                  t * (F[5] + t * (F[6] - t * F[7]))))))) /
             (F[8] +
              t * (F[9] +
                   t * (F[10] + t * (F[11] + t * (F[12] + t * (F[13] + t)))))) *
             F[14]; // correction to remove bias
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;

        fd = (G[0] +
              t * (G[1] +
                   t * (G[2] +
                        t * (G[3] +
                             t * (G[4] +
                                  t * (G[5] + t * (G[6] + t * G[7]))))))) /
             (G[8] +
              t * (G[9] +
                   t * (G[10] +
                        t * (G[11] +
                             t * (G[12] + t * (G[13] + t * (G[14] + t))))))) *
             G[15]; // correction to remove bias
    }
    else
    {
        w  = one / (x * x);
        t  = onethousandsixhundred * w;
        fd = x * x * x * x * sqrt(x) * factor *
             (H[0] +
              w * (H[1] + t * (H[2] + t * (H[3] + t * (H[4] - t * H[5])))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 4 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p4(const BS_REAL y)
{
    BS_REAL x, ex, t, s, y2;
    BS_REAL fd = 0.;

    x = -fabs(y);

    constexpr BS_REAL A[8]  = {24,
                               18247.2542465629199,
                               2120.56302902849207,
                               54.659507299984584,
                               0.121876197098273914,
                               24329.6723287505678,
                               3261.01909656925521,
                               117.071576489684022};
    constexpr BS_REAL B[15] = {
        33093.9102025608137, 19031.4783357798975, 5431.08357152245897,
        2393.94262931609398, 1268.40017257978070, 418.662172475927160,
        77.4108960876508190, 6.67374311450268357, 2645.4885670047153,
        3236.2237702166948,  2500.5977847175497,  1278.0109577275445,
        448.99020896813485,  105.86020755838874,  15.216887271751039};
    constexpr BS_REAL C[3] = {45.4575758158678040, 6.57973626739290575, 0.2};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex * (A[0] - ex * (A[1] + t * (A[2] + t * (A[3] + t * A[4]))) /
                              (A[5] + t * (A[6] + t * (A[7] + t))));
    }
    else if (x <= zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] +
                             t * (B[4] +
                                  t * (B[5] + t * (B[6] + t * B[7]))))))) /
             (B[8] +
              s * (B[9] +
                   s * (B[10] +
                        s * (B[11] +
                             s * (B[12] + s * (B[13] + s * (B[14] + s)))))));
    }
    if (y > zero)
    {
        y2 = y * y;
        fd = fd + y * (C[0] + y2 * (C[1] + y2 * C[2]));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 9/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p92(const BS_REAL x)
{
    const BS_REAL factor = 2. / 11.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[8]  = {52.3427777845535202,  37544.8291647986266,
                               4113.48125371067712,  99.866490020209337,
                               0.221035794126711913, 32460.7344732907336,
                               4028.81223496037193,  131.831032502212893};
    constexpr BS_REAL B[15] = {
        69576.0137991453185, 40945.8214269928341, 11636.1644381151941,
        4627.86250549447612, 2350.13392862571985, 772.917515779084181,
        143.282768049011730, 12.4161720203166775, 2535.8419026162019,
        3176.5677917330904,  2403.1400307314004,  1205.9132705211921,
        416.75320914741449,  97.44802030226960,   14.192502951043942};
    constexpr BS_REAL C[15] = {
        514822.53904173208,  522785.81416122782,  325144.232801546845,
        158667.557186950475, 59144.454873935605,  15870.4677825005999,
        2843.37449801975925, 272.600692278110437, 4614.934746888508,
        2755.8568288489921,  1644.4943294235255,  712.2181577760027,
        238.45583147179604,  59.96344339125020,   10.461476374827578};
    constexpr BS_REAL D[15] = {
        3.42730803487596668e7, 9.75806658972896811e7, 1.34864192616165924e8,
        1.16686466196565430e8, 6.80744250678140765e7, 2.67765233037235380e7,
        6.59736845425507999e6, 792820.676193592762,   98882.7648048876841,
        8989.62788188007021,   9880.17356654342362,   -2706.34081341857522,
        689.741654284866509,   -130.873077409052406,  16.2961910923971688};
    constexpr BS_REAL E[16] = {
        4.89775289337907328e9,  2.18308560190839931e10, 4.37929196886803725e10,
        5.12376841913497741e10, 3.77548677711700977e10, 1.75249365659518690e10,
        4.74974437725451705e9,  5.81264585063986499e8,  1.16603582458148748e6,
        904387.524377664561,    128592.404626994618,    -10482.2282242040897,
        1464.66129358886087,    -193.730497094764856,   19.0689920747890676,
        0.999999999999999823};
    constexpr BS_REAL F[16] = {
        2.2676209606423816e13, 1.24113302112340516e14, 2.9428117821474437e14,
        3.9045588584587502e14, 3.1161725148614816e14,  1.4855808500445768e14,
        3.8726601091729722e13, 4.1673615587969372e12,  2.72286381952640596e8,
        1.76997321474048891e8, 1.58695998922133285e7,  -757839.84310180434,
        61065.315196529561,    -4351.64276565798757,   174.373773100542411,
        0.999999999999999697};
    constexpr BS_REAL G[16] = {
        1.65627579413599954e13, 9.8868092265530109e13,  2.51632614382698224e14,
        3.53003554017750602e14, 2.93707971643967604e14, 1.44188847117098389e14,
        3.83794461491062424e13, 4.22247796386900843e12, 5.7640009494074862e6,
        3.82118971474499980e6,  375185.320814949919,    -21189.7277758964890,
        2225.13154291024219,    -240.494104779578652,   20.7062411414387229,
        0.999999999999999697};
    constexpr BS_REAL H[6] = {1,
                              40.7121181544936151,
                              0.256364746421960192,
                              0.000125058641885852090,
                              5.92447825215879480e-8,
                              3.38237202703194112e-10};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex * (A[0] - ex * (A[1] + t * (A[2] + t * (A[3] + t * A[4]))) /
                              (A[5] + t * (A[6] + t * (A[7] + t))));
    }
    else if (x < zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] +
                             t * (B[4] +
                                  t * (B[5] + t * (B[6] + t * B[7]))))))) /
             (B[8] +
              s * (B[9] +
                   s * (B[10] +
                        s * (B[11] +
                             s * (B[12] + s * (B[13] + s * (B[14] + s)))))));
    }
    else if (x < two)
    {
        t = half * x;
        s = one - t;

        fd = (C[0] +
              t * (C[1] +
                   t * (C[2] +
                        t * (C[3] +
                             t * (C[4] +
                                  t * (C[5] + t * (C[6] + t * C[7]))))))) /
             (C[8] +
              s * (C[9] +
                   s * (C[10] +
                        s * (C[11] +
                             s * (C[12] + s * (C[13] + s * (C[14] + s)))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);

        fd = (D[0] +
              t * (D[1] +
                   t * (D[2] +
                        t * (D[3] +
                             t * (D[4] +
                                  t * (D[5] + t * (D[6] + t * D[7]))))))) /
             (D[8] +
              t * (D[9] +
                   t * (D[10] +
                        t * (D[11] +
                             t * (D[12] + t * (D[13] + t * (D[14] - t)))))));
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;

        fd = (E[0] +
              t * (E[1] +
                   t * (E[2] +
                        t * (E[3] +
                             t * (E[4] +
                                  t * (E[5] + t * (E[6] + t * E[7]))))))) /
             (E[8] +
              t * (E[9] +
                   t * (E[10] +
                        t * (E[11] +
                             t * (E[12] + t * (E[13] + t * (E[14] - t))))))) *
             E[15]; // correction to remove bias
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;

        fd = (F[0] +
              t * (F[1] +
                   t * (F[2] +
                        t * (F[3] +
                             t * (F[4] +
                                  t * (F[5] + t * (F[6] + t * F[7]))))))) /
             (F[8] +
              t * (F[9] +
                   t * (F[10] +
                        t * (F[11] +
                             t * (F[12] + t * (F[13] + t * (F[14] + t))))))) *
             F[15]; // correction to remove bias
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;

        fd = (G[0] +
              t * (G[1] +
                   t * (G[2] +
                        t * (G[3] +
                             t * (G[4] +
                                  t * (G[5] + t * (G[6] + t * G[7]))))))) /
             (G[8] +
              t * (G[9] +
                   t * (G[10] +
                        t * (G[11] +
                             t * (G[12] + t * (G[13] + t * (G[14] - t))))))) *
             G[15]; // correction to remove bias
    }
    else
    {
        w = one / (x * x);
        t = onethousandsixhundred * w;

        fd = x * x * x * x * x * sqrt(x) * factor *
             (H[0] +
              w * (H[1] + t * (H[2] + t * (H[3] + t * (H[4] + t * H[5])))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 5 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p5(const BS_REAL y)
{
    BS_REAL x, ex, t, s, y2;
    BS_REAL fd = 0.;

    x = -fabs(y);

    constexpr BS_REAL A[8]  = {120,
                               81190.938912603315,
                               8368.4990332049831,
                               190.753830813273698,
                               0.413800735538960261,
                               43301.8340867217726,
                               4977.68099709243407,
                               148.484432990224213};
    constexpr BS_REAL B[15] = {
        159651.547840031111, 96307.2005742063042, 26923.5693307648389,
        9274.54751848368696, 4445.76333033698006, 1461.45267097859337,
        272.164427980501432, 23.6526046298891619, 2522.7839609396783,
        3244.5527999477567,  2403.0532924519756,  1176.7202478443275,
        397.7596246691212,   91.84661231161838,   13.491911254479298};
    constexpr BS_REAL C[4] = {236.532261911384425, 113.643939539669510,
                              8.22467033424113218, 0.166666666666666667};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex * (A[0] - ex * (A[1] + t * (A[2] + t * (A[3] + t * A[4]))) /
                              (A[5] + t * (A[6] + t * (A[7] + t))));
    }
    else if (x <= zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] +
                             t * (B[4] +
                                  t * (B[5] + t * (B[6] + t * B[7]))))))) /
             (B[8] +
              s * (B[9] +
                   s * (B[10] +
                        s * (B[11] +
                             s * (B[12] + s * (B[13] + s * (B[14] + s)))))));
    }
    if (y > zero)
    {
        y2 = y * y;

        fd = -fd + C[0] + y2 * (C[1] + y2 * (C[2] + y2 * C[3]));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 11/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p112(const BS_REAL x)
{
    const BS_REAL factor = 2. / 13.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[8]  = {287.885277815044361, 183706.438585548927,
                               17781.2198813711389, 379.468701417596835,
                               0.79823396490661890, 57756.370490854840,
                               6150.6214356632731,  167.282750039076672};
    constexpr BS_REAL B[15] = {
        399729.728975269086, 245048.290430557943, 64937.4257386015642,
        17956.2886271385592, 8102.86054746075370, 2739.56619011305740,
        522.197272134368319, 45.9656786040236010, 2594.1469579416305,
        3458.1847294224010,  2522.2748339106920,  1201.9024264852580,
        394.94834682669296,  89.38046041124415,   13.110503750759739};
    constexpr BS_REAL C[15] = {
        2.40644550206125935e6, 2.04903536564989627e6, 1.06804765545666824e6,
        478219.980682203232,   171969.059741023943,   43927.0876916260199,
        7256.0914773211089,    626.34457422853533,    3109.1225672233219,
        2636.2757909353112,    1626.1798939585643,    743.1451899485297,
        255.41213934754826,    64.38943570079975,     10.994511861529934};
    constexpr BS_REAL D[17] = {
        7.70852968624936390e8, 2.17991233625412353e9, 3.06662111607913460e9,
        2.78373598102826860e9, 1.77557231061500545e9, 8.12988256861068453e8,
        2.60886497846886753e8, 5.37757509332685362e7, 5.49029147083874521e6,
        384975.653810726401,   -10859.8352547471390,  47234.6496298591730,
        -15689.8880335197160,  4580.68161072037214,   -1041.20834072320464,
        172.979919201040627,   -18.7593019718622089};
    constexpr BS_REAL E[16] = {
        1.71681959469725307e10, 7.11557910478070659e10, 1.35103963154256158e11,
        1.52010553996923224e11, 1.09288906897533869e11, 5.01817283608162539e10,
        1.36404590384355901e10, 1.69961340800732270e9,  600169.964027912926,
        63989.8050690151317,    3306.49827047230707,    -2248.31637427985521,
        649.075943625483938,    -126.914632653258205,   16.0289311039308149,
        0.999999999999999735};
    constexpr BS_REAL F[16] = {
        1.05700178515928121e12, 6.12117465680051512e12, 1.55582880442596747e13,
        2.24736526624247277e13, 1.99069164118763497e13, 1.08056362404911750e13,
        3.32622572257609455e12, 4.47682895372084249e11, 1.27888301811042482e6,
        318603.811469994946,    -33834.0645458632405,   6000.96391500182524,
        -1070.68701259288375,   161.612080255393227,    -17.5331213409665917,
        0.999999999999999653};
    constexpr BS_REAL G[16] = {
        1.06038260567473087e14, 6.8238792169646891e14,  1.88917697801412952e15,
        2.91401892625408501e15, 2.70150682400302739e15, 1.50286942831009197e15,
        4.63446459437813772e14, 6.0885766731916815e13,  2.08723140369481691e6,
        445823.690247689681,    -42239.1125174643989,   6891.8137509933955,
        -1158.79785029794968,   168.017431473136672,    -17.7743748531881933,
        0.999999999999999637};
    constexpr BS_REAL H[6] = {1,
                              58.8063928898240827,
                              0.666548340699108862,
                              0.00162576228417182612,
                              -2.56764926438640722e-7,
                              6.18964665583859548e-10};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex * (A[0] - ex * (A[1] + t * (A[2] + t * (A[3] + t * A[4]))) /
                              (A[5] + t * (A[6] + t * (A[7] + t))));
    }
    else if (x < zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] +
                             t * (B[4] +
                                  t * (B[5] + t * (B[6] + t * B[7]))))))) /
             (B[8] +
              s * (B[9] +
                   s * (B[10] +
                        s * (B[11] +
                             s * (B[12] + s * (B[13] + s * (B[14] + s)))))));
    }
    else if (x < two)
    {
        t = half * x;
        s = one - t;

        fd = (C[0] +
              t * (C[1] +
                   t * (C[2] +
                        t * (C[3] +
                             t * (C[4] +
                                  t * (C[5] + t * (C[6] + t * C[7]))))))) /
             (C[8] +
              s * (C[9] +
                   s * (C[10] +
                        s * (C[11] +
                             s * (C[12] + s * (C[13] + s * (C[14] + s)))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);

        fd =
            (D[0] +
             t * (D[1] +
                  t * (D[2] +
                       t * (D[3] +
                            t * (D[4] +
                                 t * (D[5] +
                                      t * (D[6] + t * (D[7] + t * D[8])))))))) /
            (D[9] +
             t * (D[10] +
                  t * (D[11] +
                       t * (D[12] +
                            t * (D[13] +
                                 t * (D[14] +
                                      t * (D[15] + t * (D[16] + t))))))));
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;

        fd = (E[0] +
              t * (E[1] +
                   t * (E[2] +
                        t * (E[3] +
                             t * (E[4] +
                                  t * (E[5] + t * (E[6] + t * E[7]))))))) /
             (E[8] +
              t * (E[9] +
                   t * (E[10] +
                        t * (E[11] +
                             t * (E[12] + t * (E[13] + t * (E[14] - t))))))) *
             E[15];
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;

        fd = (F[0] +
              t * (F[1] +
                   t * (F[2] +
                        t * (F[3] +
                             t * (F[4] +
                                  t * (F[5] + t * (F[6] + t * F[7]))))))) /
             (F[8] +
              t * (F[9] +
                   t * (F[10] +
                        t * (F[11] +
                             t * (F[12] + t * (F[13] + t * (F[14] + t))))))) *
             F[15];
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;

        fd = (G[0] +
              t * (G[1] +
                   t * (G[2] +
                        t * (G[3] +
                             t * (G[4] +
                                  t * (G[5] + t * (G[6] + t * G[7]))))))) /
             (G[8] +
              t * (G[9] +
                   t * (G[10] +
                        t * (G[11] +
                             t * (G[12] + t * (G[13] + t * (G[14] + t))))))) *
             G[15];
    }
    else
    {
        w = one / (x * x);
        t = onethousandsixhundred * w;

        fd = x * x * x * x * x * x * sqrt(x) * factor *
             (H[0] +
              w * (H[1] + t * (H[2] + t * (H[3] + t * (H[4] - t * H[5])))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 6 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p6(const BS_REAL y)
{
    BS_REAL x, ex, t, s, y2;
    BS_REAL fd = 0.;

    x = -fabs(y);

    constexpr BS_REAL A[8]  = {720,
                               433290.557356514403,
                               39323.1678283287030,
                               783.71684376947655,
                               1.58412947146158337,
                               77029.432418935897,
                               7600.9245809611507,
                               188.511069473956679};
    constexpr BS_REAL B[15] = {
        1.06809887479312876e6, 651074.246191348755, 152197.899924352192,
        27269.0203707062592,   11937.9088600476726, 4659.76125900467198,
        964.100791156809939,   88.5841245838029230, 2681.3731718905701,
        3764.9048490469408,    2739.9504946219358,  1276.0933863294022,
        406.93880737411049,    89.70143752313466,   12.997281214703279};
    constexpr BS_REAL C[4] = {1419.19357146830655, 227.287879079339020,
                              9.86960440108935862, 0.142857142857142857};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex * (A[0] - ex * (A[1] + t * (A[2] + t * (A[3] + t * A[4]))) /
                              (A[5] + t * (A[6] + t * (A[7] + t))));
    }
    else if (x <= zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] +
                             t * (B[4] +
                                  t * (B[5] + t * (B[6] + t * B[7]))))))) /
             (B[8] +
              s * (B[9] +
                   s * (B[10] +
                        s * (B[11] +
                             s * (B[12] + s * (B[13] + s * (B[14] + s)))))));
    }
    if (y > zero)
    {
        y2 = y * y;

        fd = fd + y * (C[0] + y2 * (C[1] + y2 * (C[2] + y2 * C[3])));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 13/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p132(const BS_REAL x)
{
    const BS_REAL factor = 2. / 15.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[7]  = {1871.25430579778835,      10.2714824358192133,
                               0.0648393897767320738,    0.000971410002224865156,
                               0.0000232470542763426993, 7.41058598210549775e-7,
                               3.37424008840868627e-8};
    constexpr BS_REAL B[15] = {
        2.94504307654357499e6, 1.72142318914178000e6, 302664.670131641659,
        -4392.91160671565775,  1886.14265293074511,   5292.33536016064607,
        1532.39748567462704,   160.488909859013319,   2672.0262186872157,
        4050.0534640335863,    2998.5673483066873,    1381.1929813352472,
        429.69401904358381,    92.16501861997658,     13.091741877025778};
    constexpr BS_REAL C[15] = {
        1.45543624525939873e7, 1.09901378809591560e7, 5.0278309600981774e6,
        2.17920053476107569e6, 801203.84933034017,    204535.562039585858,
        32251.6713113453623,   2520.49966644549104,   2525.1970865635211,
        2519.0633694107075,    1640.0689226392724,    779.5069098915619,
        273.74080946976357,    69.07956623994342,     11.541634429054885};
    constexpr BS_REAL D[17] = {
        1.93509493140212933e9, 4.73031979570464235e9, 5.87850060979944364e9,
        4.82135265081185025e9, 2.83819099819062720e9, 1.22325945046779759e9,
        3.76487745749965788e8, 7.58303164598034022e7, 7.71973955404921694e6,
        101717.15224358499,    25406.290690742596,    11612.352340539603,
        4174.0221825716411,    1281.7558184421595,    331.85386698947176,
        69.634528227703809,    10.868655326785347};
    constexpr BS_REAL E[18] = {
        5.2085062280477356e11,  2.2833424052547072e12, 4.6942573239060036e12,
        5.8958834991587034e12,  4.9360061465395771e12, 2.8153080727145281e12,
        1.06730054566227826e12, 2.4600629203104023e11, 2.6466259504577480e10,
        2.42061644583385800e6,  153096.878206723536,   38844.6431417054012,
        -16885.0152281400560,   4865.05644899569952,   -1067.44742924509046,
        173.008767319757249,    -18.6099093504480426,  0.999999999999999770};
    constexpr BS_REAL F[18] = {
        7.7494086475291441e13,  5.0061227151520991e14,  1.45146838868293917e15,
        2.46379854838342831e15, 2.67523099025981612e15, 1.90126580903196655e15,
        8.6325468139568610e14,  2.28901218799706877e14, 2.71416555283189345e13,
        8.9593891719102305e6,   2.23010123362691961e6,  -239561.734722929070,
        43938.8978835860601,    -8408.3951285804076,    1448.44898008864648,
        -201.658017368512975,   19.7129633131788616,    0.999999999999999586};
    constexpr BS_REAL G[18] = {
        1.22298919425291664e16, 8.9726625177108898e16, 2.89429751278928300e17,
        5.3577901105140076e17,  6.2206348573306568e17, 4.63426239396617756e17,
        2.16065376561259960e17, 5.7540345485220279e16, 6.6835199866465357e15,
        1.31871866498468464e7,  2.83885323472743888e6, -275298.111415499841,
        47078.0598194149218,    -8612.3789225919783,   1446.29041011108831,
        -199.315359842233309,   19.5132779151984627,   0.999999999999999580};
    constexpr BS_REAL H[6] = {1,
                              80.1905357588510589,
                              1.42831787292433216,
                              0.00812881145350198890,
                              3.85165680520411590e-6,
                              1.83598087001386478e-9};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        s = one - t;
        fd =
            ex * (A[0] -
                  ex * (A[1] +
                        s * (A[2] +
                             s * (A[3] + s * (A[4] + s * (A[5] + s * A[6]))))));
    }
    else if (x < zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] +
                             t * (B[4] +
                                  t * (B[5] + t * (B[6] + t * B[7]))))))) /
             (B[8] +
              s * (B[9] +
                   s * (B[10] +
                        s * (B[11] +
                             s * (B[12] + s * (B[13] + s * (B[14] + s)))))));
    }
    else if (x < two)
    {
        t = half * x;
        s = one - t;

        fd = (C[0] +
              t * (C[1] +
                   t * (C[2] +
                        t * (C[3] +
                             t * (C[4] +
                                  t * (C[5] + t * (C[6] + t * C[7]))))))) /
             (C[8] +
              s * (C[9] +
                   s * (C[10] +
                        s * (C[11] +
                             s * (C[12] + s * (C[13] + s * (C[14] + s)))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);
        s = one - t;

        fd =
            (D[0] +
             t * (D[1] +
                  t * (D[2] +
                       t * (D[3] +
                            t * (D[4] +
                                 t * (D[5] +
                                      t * (D[6] + t * (D[7] + t * D[8])))))))) /
            (D[9] +
             s * (D[10] +
                  s * (D[11] +
                       s * (D[12] +
                            s * (D[13] +
                                 s * (D[14] +
                                      s * (D[15] + s * (D[16] + s))))))));
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;

        fd =
            (E[0] +
             t * (E[1] +
                  t * (E[2] +
                       t * (E[3] +
                            t * (E[4] +
                                 t * (E[5] +
                                      t * (E[6] + t * (E[7] + t * E[8])))))))) /
            (E[9] +
             t * (E[10] +
                  t * (E[11] +
                       t * (E[12] +
                            t * (E[13] +
                                 t * (E[14] +
                                      t * (E[15] + t * (E[16] + t)))))))) *
            E[17]; // correction to remove bias
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;

        fd =
            (F[0] +
             t * (F[1] +
                  t * (F[2] +
                       t * (F[3] +
                            t * (F[4] +
                                 t * (F[5] +
                                      t * (F[6] + t * (F[7] + t * F[8])))))))) /
            (F[9] +
             t * (F[10] +
                  t * (F[11] +
                       t * (F[12] +
                            t * (F[13] +
                                 t * (F[14] +
                                      t * (F[15] + t * (F[16] - t)))))))) *
            F[17]; // correction to remove bias
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;

        fd =
            (G[0] +
             t * (G[1] +
                  t * (G[2] +
                       t * (G[3] +
                            t * (G[4] +
                                 t * (G[5] +
                                      t * (G[6] + t * (G[7] + t * G[8])))))))) /
            (G[9] +
             t * (G[10] +
                  t * (G[11] +
                       t * (G[12] +
                            t * (G[13] +
                                 t * (G[14] +
                                      t * (G[15] + t * (G[16] - t)))))))) *
            G[17]; // correction to remove bias
    }
    else
    {
        w = one / (x * x);
        t = onethousandsixhundred * w;

        fd = x * x * x * x * x * x * x * sqrt(x) * factor *
             (H[0] +
              w * (H[1] + t * (H[2] + t * (H[3] + t * (H[4] + t * H[5])))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 7 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p7(const BS_REAL y)
{
    BS_REAL x, ex, t, s, y2;
    BS_REAL fd = 0.;

    x = -fabs(y);

    constexpr BS_REAL A[7]  = {5040,
                               19.5849162780581217,
                               0.101236264600248646,
                               0.00131827373096852460,
                               0.0000283211236200358235,
                               8.27819847443616991e-7,
                               3.46664379401620231e-8};
    constexpr BS_REAL B[15] = {
        8.05214238081846197e6, 4.23358827341165564e6, 262627.912342619479,
        -323968.614686001584,  -98175.8116505823446,  -7562.84276896647246,
        1252.57372771721279,   240.247242901336449,   2413.8835947192718,
        4139.8252477104615,    3204.5076382872745,    1486.4720567970787,
        456.50435752749115,    95.80622002968487,     13.317952574698873};
    constexpr BS_REAL C[5] = {10042.0286586746908, 4967.17750013907292,
                              397.753788388843285, 11.5145384679375851, 0.125};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;
        s  = one - t;

        fd =
            ex * (A[0] -
                  ex * (A[1] +
                        s * (A[2] +
                             s * (A[3] + s * (A[4] + s * (A[5] + s * A[6]))))));
    }
    else if (x <= zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] +
                             t * (B[4] +
                                  t * (B[5] + t * (B[6] + t * B[7]))))))) /
             (B[8] +
              s * (B[9] +
                   s * (B[10] +
                        s * (B[11] +
                             s * (B[12] + s * (B[13] + s * (B[14] + s)))))));
    }
    if (y > zero)
    {
        y2 = y * y;

        fd = -fd + C[0] + y2 * (C[1] + y2 * (C[2] + y2 * (C[3] + y2 * C[4])));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 15/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p152(const BS_REAL x)
{
    const BS_REAL factor = 2. / 17.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[6]  = {14034.4072934834126, 171704.258234536602,
                               6049.5894628442878,  18.1763510883857458,
                               4429.36992795387362, 175.155804618692908};
    constexpr BS_REAL B[15] = {
        2.06790092766575034e7,  8.26618273142399596e6, -2.08115753677480168e6,
        -2.17301177649999531e6, -641552.116750411401,  -91807.9577537648588,
        -5144.12042361256270,   130.857943513312398,   1711.2966803109211,
        3804.4823677048040,     3232.4717528599742,    1550.7278916841484,
        478.50436506151174,     99.39616883636118,     13.585756243903013};
    constexpr BS_REAL C[15] = {
        1.14403769376860621e8, 8.4972567501965974e7, 3.81753609315202610e7,
        1.65881819811148351e7, 6.1705051839981204e6, 1.58079386467916017e6,
        246507.752966115175,   18602.2895319461278,  2572.2942565752562,
        2659.5730458609857,    1732.1057976186518,   829.4655494427400,
        293.34263153462293,    73.73539965269387,    12.058936149293461};
    constexpr BS_REAL D[17] = {
        6.26769725573910255e9,  1.32711699403202671e10, 1.45560142717062095e10,
        1.07860795716164403e10, 5.86969120046064229e9,  2.38796511791514625e9,
        7.07324240427174132e8,  1.39791587840460999e8,  1.42694865767192645e7,
        30361.411562018280,     17707.813587539272,     8523.204548822773,
        3390.1996842891790,     1127.1366883800605,     309.69425190845478,
        67.751604120864452,     10.834010773474702};
    constexpr BS_REAL E[18] = {
        5.69341524462926732e11, 2.36489781301030086e12, 4.67723648136547218e12,
        5.73677938363737470e12, 4.75966348844286492e12, 2.72977815388509051e12,
        1.05619383963810221e12, 2.52484808268258635e11, 2.87129101561538399e10,
        233403.35536095456,     60565.447915936497,     18770.529358157892,
        5494.0347508271098,     1458.2356321178129,     341.00432855822186,
        67.360611688552025,     10.355894270482783,     0.999999999999999687};
    constexpr BS_REAL F[18] = {
        6.79824678667324436e13, 4.47889271605157194e14, 1.33539302658610566e15,
        2.35081455707368313e15, 2.67105097393954020e15, 2.00614187037046236e15,
        9.73640078908216993e14, 2.79804648854039015e14, 3.66107518187916541e13,
        585167.34458490542,     99101.77544747319,      23659.896891961921,
        5847.0608231871316,     1394.1981365390388,     307.36846956466131,
        59.734384735038341,     9.438593558054937,      0.999999999999999620};
    constexpr BS_REAL G[18] = {
        1.99221037305471970e16, 1.55301826780028291e17, 5.3562790972794281e17,
        1.06787028331863552e18, 1.34677794609411928e18, 1.10118289362700600e18,
        5.7076444557408842e17,  1.71793041083229002e17, 2.30600690809359773e16,
        959281.46439522639,     146009.54484661068,     31788.065669901676,
        7257.6052635887384,     1618.0190118885021,     337.33224463190506,
        62.698857808593658,     9.5893450572160768,     0.999999999999999649};
    constexpr BS_REAL H[6] = {1,
                              104.864546761574387,
                              2.69793375997248149,
                              0.0276379588953803004,
                              0.0000654786006189425666,
                              1.03209499898071117e-8};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex * (A[0] - ex * (A[1] + t * (A[2] + t * A[3])) /
                              (A[4] + t * (A[5] + t)));
    }
    else if (x < zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] +
                             t * (B[4] +
                                  t * (B[5] + t * (B[6] + t * B[7]))))))) /
             (B[8] +
              s * (B[9] +
                   s * (B[10] +
                        s * (B[11] +
                             s * (B[12] + s * (B[13] + s * (B[14] + s)))))));
    }
    else if (x < two)
    {
        t = half * x;
        s = one - t;

        fd = (C[0] +
              t * (C[1] +
                   t * (C[2] +
                        t * (C[3] +
                             t * (C[4] +
                                  t * (C[5] + t * (C[6] + t * C[7]))))))) /
             (C[8] +
              s * (C[9] +
                   s * (C[10] +
                        s * (C[11] +
                             s * (C[12] + s * (C[13] + s * (C[14] + s)))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);
        s = one - t;

        fd =
            (D[0] +
             t * (D[1] +
                  t * (D[2] +
                       t * (D[3] +
                            t * (D[4] +
                                 t * (D[5] +
                                      t * (D[6] + t * (D[7] + t * D[8])))))))) /
            (D[9] +
             s * (D[10] +
                  s * (D[11] +
                       s * (D[12] +
                            s * (D[13] +
                                 s * (D[14] +
                                      s * (D[15] + s * (D[16] + s))))))));
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;
        s = one - t;

        fd =
            (E[0] +
             t * (E[1] +
                  t * (E[2] +
                       t * (E[3] +
                            t * (E[4] +
                                 t * (E[5] +
                                      t * (E[6] + t * (E[7] + t * E[8])))))))) /
            (E[9] +
             s * (E[10] +
                  s * (E[11] +
                       s * (E[12] +
                            s * (E[13] +
                                 s * (E[14] +
                                      s * (E[15] + s * (E[16] + s)))))))) *
            E[17]; // correction to remove bias
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;
        s = one - t;

        fd =
            (F[0] +
             t * (F[1] +
                  t * (F[2] +
                       t * (F[3] +
                            t * (F[4] +
                                 t * (F[5] +
                                      t * (F[6] + t * (F[7] + t * F[8])))))))) /
            (F[9] +
             s * (F[10] +
                  s * (F[11] +
                       s * (F[12] +
                            s * (F[13] +
                                 s * (F[14] +
                                      s * (F[15] + s * (F[16] + s)))))))) *
            F[17]; // correction to remove bias
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;
        s = one - t;

        fd =
            (G[0] +
             t * (G[1] +
                  t * (G[2] +
                       t * (G[3] +
                            t * (G[4] +
                                 t * (G[5] +
                                      t * (G[6] + t * (G[7] + t * G[8])))))))) /
            (G[9] +
             s * (G[10] +
                  s * (G[11] +
                       s * (G[12] +
                            s * (G[13] +
                                 s * (G[14] +
                                      s * (G[15] + s * (G[16] + s)))))))) *
            G[17]; // correction to remove bias
    }
    else
    {
        w = one / (x * x);
        t = onethousandsixhundred * w;

        fd = x * x * x * x * x * x * x * x * sqrt(x) * factor *
             (H[0] +
              w * (H[1] + t * (H[2] + t * (H[3] + t * (H[4] - t * H[5])))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 8 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p8(const BS_REAL y)
{
    BS_REAL x, ex, t, s, y2;
    BS_REAL fd = 0.;

    x = -fabs(y);

    constexpr BS_REAL A[6]  = {40320,
                               438381.668209835602,
                               14167.5566380822417,
                               39.1239843501671193,
                               5566.7513423478282,
                               199.502568378336498};
    constexpr BS_REAL B[14] = {
        2.14460848478648315e8, 1.66513461191585055e8, 5.17491656951190442e7,
        6.60393533940475286e6, -301119.943009454383,  -202842.867158865294,
        20669.3789063593406,   10903.784789014007,    14623.305826998018,
        9179.893497142397,     3522.2940212081378,    905.1498999681050,
        159.07578981158691,    18.015010064955669};
    constexpr BS_REAL C[5] = {80336.2292693975266, 13245.8066670375278,
                              636.406061422149257, 13.1594725347858115,
                              0.111111111111111111};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;

        fd = ex * (A[0] - ex * (A[1] + t * (A[2] + t * A[3])) /
                              (A[4] + t * (A[5] + t)));
    }
    else if (x <= zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] + t * (B[4] + t * (B[5] - t * B[6])))))) /
             (B[7] +
              s * (B[8] +
                   s * (B[9] +
                        s * (B[10] +
                             s * (B[11] + s * (B[12] + s * (B[13] + s)))))));
    }
    if (y > zero)
    {
        y2 = y * y;

        fd = fd +
             y * (C[0] + y2 * (C[1] + y2 * (C[2] + y2 * (C[3] + y2 * C[4]))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 17/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p172(const BS_REAL x)
{
    const BS_REAL factor = 2. / 19.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[6]  = {119292.461994609007,      164.281538211493590,
                               0.465418918711191171,     0.00397356161994318188,
                               0.0000615968461595874662, 1.50191715167058393e-6};
    constexpr BS_REAL B[13] = {
        1.43233835112073902e8, 1.39644553704144873e8, 6.13766014358067700e7,
        1.56246791599180488e7, 2.43295219093432312e6, 218176.915876999296,
        8461.46979813306275,   3043.1981753546660,    3390.9687009373892,
        1761.0786885568746,    550.68360665579008,    112.08528557398394,
        14.640373761706403};
    constexpr BS_REAL C[15] = {
        1.21473360069923149e9, 1.01176690309922349e9, 4.95296056344158641e8,
        2.06442737728860908e8, 7.0192333337028334e7,  1.68572449951275985e7,
        2.53144559406161612e6, 186718.775025914870,   3455.142635040908,
        3415.073282258220,     2023.9784836523565,    903.6597935624893,
        308.93682588995655,    76.49008105409234,     12.338920623517048};
    constexpr BS_REAL D[17] = {
        2.67440099736518804e10, 4.91370104189578919e10, 4.73732500104890339e10,
        3.16238224173646457e10, 1.59135818804699047e10, 6.12061461496369111e9,
        1.74615616554871192e9,  3.38446352294444219e8,  3.46053509281020592e7,
        11121.047017722056,     9829.804681621007,      5748.503365098324,
        2606.1477102230902,     953.9183383047316,      281.64461270648236,
        64.892044676666959,     10.720667663673541};
    constexpr BS_REAL E[18] = {
        1.21689593180741239e12, 4.79662223812667505e12, 9.12838410654446290e12,
        1.09342887556403655e13, 8.99572172197880880e12, 5.19592600765981107e12,
        2.05751238857523311e12, 5.12310527144160372e11, 6.19902241659841813e10,
        39255.392055313674,     22199.471524635690,     9415.210433779872,
        3354.6940213500329,     1027.9861625160893,     269.15608801102034,
        58.326091875891505,     9.690839898241749,      0.999999999999999611};
    constexpr BS_REAL F[18] = {
        8.36136698862497927e13, 5.63471803594835697e14, 1.73344511753855139e15,
        3.17713365223451161e15, 3.79488820022826435e15, 3.02846836744419910e15,
        1.58129984938571763e15, 4.96500350365744464e14, 7.24731492178964020e13,
        43633.188864890011,     21099.956516665469,     8040.798889514530,
        2678.4837772160305,     794.3651890403060,      208.01106751201808,
        46.655407061024955,     8.3567416627775934,     0.999999999999999524};
    constexpr BS_REAL G[18] = {
        2.85846778936881832e16, 2.36061789976170730e17, 8.67582358035872224e17,
        1.85563671344211925e18, 2.53068942615802160e18, 2.25901769990104319e18,
        1.29359182644194722e18, 4.36848670280351452e17, 6.72189933519645798e16,
        50522.192743560920,     23134.306793549278,     8428.128974587599,
        2708.9200442813342,     782.3736655720094,      201.43688280680891,
        44.899007032888684,     8.0956175513436493,     0.999999999999999478};
    constexpr BS_REAL H[6] = {1,
                              132.828425897994463,
                              4.66006740357572083,
                              0.0750173171258174240,
                              0.000414697262104018271,
                              1.95167395079296505e-7};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;
        s  = one - t;

        fd = ex *
             (A[0] -
              ex * (A[1] + s * (A[2] + s * (A[3] + s * (A[4] + s * A[5])))));
    }
    else if (x < zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] + t * (B[4] + t * (B[5] + t * B[6])))))) /
             (B[7] +
              s * (B[8] +
                   s * (B[9] + s * (B[10] + s * (B[11] + s * (B[12] + s))))));
    }
    else if (x < two)
    {
        t = half * x;
        s = one - t;

        fd = (C[0] +
              t * (C[1] +
                   t * (C[2] +
                        t * (C[3] +
                             t * (C[4] +
                                  t * (C[5] + t * (C[6] + t * C[7]))))))) /
             (C[8] +
              s * (C[9] +
                   s * (C[10] +
                        s * (C[11] +
                             s * (C[12] + s * (C[13] + s * (C[14] + s)))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);
        s = one - t;

        fd =
            (D[0] +
             t * (D[1] +
                  t * (D[2] +
                       t * (D[3] +
                            t * (D[4] +
                                 t * (D[5] +
                                      t * (D[6] + t * (D[7] + t * D[8])))))))) /
            (D[9] +
             s * (D[10] +
                  s * (D[11] +
                       s * (D[12] +
                            s * (D[13] +
                                 s * (D[14] +
                                      s * (D[15] + s * (D[16] + s))))))));
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;
        s = one - t;

        fd =
            (E[0] +
             t * (E[1] +
                  t * (E[2] +
                       t * (E[3] +
                            t * (E[4] +
                                 t * (E[5] +
                                      t * (E[6] + t * (E[7] + t * E[8])))))))) /
            (E[9] +
             s * (E[10] +
                  s * (E[11] +
                       s * (E[12] +
                            s * (E[13] +
                                 s * (E[14] +
                                      s * (E[15] + s * (E[16] + s)))))))) *
            E[17];
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;
        s = one - t;

        fd =
            (F[0] +
             t * (F[1] +
                  t * (F[2] +
                       t * (F[3] +
                            t * (F[4] +
                                 t * (F[5] +
                                      t * (F[6] + t * (F[7] + t * F[8])))))))) /
            (F[9] +
             s * (F[10] +
                  s * (F[11] +
                       s * (F[12] +
                            s * (F[13] +
                                 s * (F[14] +
                                      s * (F[15] + s * (F[16] + s)))))))) *
            F[17];
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;
        s = one - t;

        fd =
            (G[0] +
             t * (G[1] +
                  t * (G[2] +
                       t * (G[3] +
                            t * (G[4] +
                                 t * (G[5] +
                                      t * (G[6] + t * (G[7] + t * G[8])))))))) /
            (G[9] +
             s * (G[10] +
                  s * (G[11] +
                       s * (G[12] +
                            s * (G[13] +
                                 s * (G[14] +
                                      s * (G[15] + s * (G[16] + s)))))))) *
            G[17];
    }
    else
    {
        w = one / (x * x);
        t = onethousandsixhundred * w;

        fd = x * x * x * x * x * x * x * x * x * sqrt(x) * factor *
             (H[0] +
              w * (H[1] + t * (H[2] + t * (H[3] + t * (H[4] + t * H[5])))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 9 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p9(const BS_REAL y)
{
    BS_REAL x, ex, t, s, y2;
    BS_REAL fd = 0.;

    x = -fabs(y);

    constexpr BS_REAL A[6]  = {362880,
                               3.11619728815890900e6,
                               84459.263086739990,
                               194.729051789880448,
                               8793.5020477151129,
                               258.970670923338654};
    constexpr BS_REAL B[13] = {
        5.34284039564434465e8, 5.33773366525970850e8, 2.42092992302655038e8,
        6.42814979745526752e7, 1.06272445991848042e7, 1.04623859092301352e6,
        48187.9003961411735,   3823.5401164205447,    4167.5033491759870,
        2107.9832247353650,    639.46059745678970,    125.53444672868253,
        15.642276104269218};
    constexpr BS_REAL C[6] = {725062.913034521571, 361513.031712288870,
                              29803.0650008344375, 954.609092133223885,
                              14.8044066016340379, 0.1};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;
        fd = ex * (A[0] - ex * (A[1] + t * (A[2] + t * A[3])) /
                              (A[4] + t * (A[5] + t)));
    }
    else if (x <= zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] + t * (B[4] + t * (B[5] + t * B[6])))))) /
             (B[7] +
              s * (B[8] +
                   s * (B[9] + s * (B[10] + s * (B[11] + s * (B[12] + s))))));
    }
    if (y > zero)
    {
        y2 = y * y;

        fd = -fd + C[0] +
             y2 * (C[1] + y2 * (C[2] + y2 * (C[3] + y2 * (C[4] + y2 * C[5]))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 19/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p192(const BS_REAL x)
{
    const BS_REAL factor = 2. / 21.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[6]  = {1.13327838894878557e6,   781.077395263087427,
                               1.48017751152694790,     0.00952643972760576450,
                               0.000118887473010527875, 2.40377778104657177e-6};
    constexpr BS_REAL B[13] = {
        1.98943681340803829e9, 2.02504701261695886e9, 9.40112997876683365e8,
        2.57205491041485272e8, 4.42589158347859287e7, 4.61202797242277802e6,
        232199.034580285755,   4645.3648503361484,    4976.8851201391641,
        2464.6700288281232,    729.19053696862632,    138.78343942961328,
        16.586263194372411};
    constexpr BS_REAL C[15] = {
        2.0904667585793652e10, 1.9423314417935927e10, 9.4466078397530459e9,
        3.3380904172376415e9,  9.2799283366292994e8,  1.9279930343506925e8,
        2.6584305328117208e7,  1.8721659680393183e6,  6510.453486834102,
        6583.316420184059,     3538.873468101783,     1337.9108450581962,
        388.7430265102922,     85.60461680432659,     12.861326654789701};
    constexpr BS_REAL D[17] = {
        1.46863114809614622e11, 2.34916998394298104e11, 1.98120796400802870e11,
        1.18887551091630769e11, 5.55928633602927965e10, 2.03490255735994097e10,
        5.60239395267172792e9,  1.05980329750305191e9,  1.07474670751003016e8,
        4912.718774356158,      5658.595947000664,      3911.882548894156,
        2003.2886387850459,     805.2407096838510,      255.52601531171211,
        62.070351460815245,     10.603837914443282};
    constexpr BS_REAL E[18] = {
        4.12413857105017004e12, 1.53777683873017139e13, 2.80178767823520436e13,
        3.25874492182871964e13, 2.64367160554342534e13, 1.53023378021225309e13,
        6.17571294140739606e12, 1.59629111145863694e12, 2.05057950138673222e11,
        9870.066781117703,      8442.670237353289,      4641.969549915204,
        1995.0747902033748,     707.3070308400987,      208.47240484486439,
        49.869236685972587,     9.0113347734666514,     0.999999999999999624};
    constexpr BS_REAL F[20] = {
        4.83390017425457148e15, 3.49375860214157535e16, 1.17174577652359819e17,
        2.38956570804562992e17, 3.26219845844040472e17, 3.09099843298003147e17,
        2.03426110846076766e17, 8.98530419590583886e16, 2.42666818156258567e16,
        3.07539397786247273e15, 207639.00554015383,     101594.67784669785,
        39318.472360728631,     13385.302295791795,     4099.2485074048775,
        1128.0104124272086,     274.41287517597127,     56.961169613859215,
        9.354997067031336,      0.999999999999999475};
    constexpr BS_REAL G[20] = {
        2.81147322924282096e18, 2.55872925207166842e19, 1.05126940075992923e20,
        2.56086191531691937e20, 4.08006931722902169e20, 4.41541112498720876e20,
        3.25219940350974164e20, 1.57672402355387744e20, 4.58547965512713670e19,
        6.13678661161381742e18, 253648.48647211156,     116226.83914995012,
        42576.299330336023,     13857.868987762036,     4097.2448088308021,
        1099.1091735270184,     263.31471897218326,     54.436947056837050,
        9.0281927555123252,     0.999999999999999480};
    constexpr BS_REAL H[7] = {1,
                              164.082173168110588,
                              7.52780119040598426,
                              0.175040406459615661,
                              0.00174172934550196671,
                              4.09368331231164635e-6,
                              6.44677689509362889e-10};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;
        s  = one - t;

        fd = ex *
             (A[0] -
              ex * (A[1] + s * (A[2] + s * (A[3] + s * (A[4] + s * A[5])))));
    }
    else if (x < zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] + t * (B[4] + t * (B[5] + t * B[6])))))) /
             (B[7] +
              s * (B[8] +
                   s * (B[9] + s * (B[10] + s * (B[11] + s * (B[12] + s))))));
    }
    else if (x < two)
    {
        t = half * x;
        s = one - t;

        fd = (C[0] +
              t * (C[1] +
                   t * (C[2] +
                        t * (C[3] +
                             t * (C[4] +
                                  t * (C[5] + t * (C[6] + t * C[7]))))))) /
             (C[8] +
              s * (C[9] +
                   s * (C[10] +
                        s * (C[11] +
                             s * (C[12] + s * (C[13] + s * (C[14] + s)))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);
        s = one - t;

        fd =
            (D[0] +
             t * (D[1] +
                  t * (D[2] +
                       t * (D[3] +
                            t * (D[4] +
                                 t * (D[5] +
                                      t * (D[6] + t * (D[7] + t * D[8])))))))) /
            (D[9] +
             s * (D[10] +
                  s * (D[11] +
                       s * (D[12] +
                            s * (D[13] +
                                 s * (D[14] +
                                      s * (D[15] + s * (D[16] + s))))))));
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;
        s = one - t;

        fd =
            (E[0] +
             t * (E[1] +
                  t * (E[2] +
                       t * (E[3] +
                            t * (E[4] +
                                 t * (E[5] +
                                      t * (E[6] + t * (E[7] + t * E[8])))))))) /
            (E[9] +
             s * (E[10] +
                  s * (E[11] +
                       s * (E[12] +
                            s * (E[13] +
                                 s * (E[14] +
                                      s * (E[15] + s * (E[16] + s)))))))) *
            E[17]; // correction to remove bias
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;
        s = one - t;

        fd = (F[0] +
              t * (F[1] +
                   t * (F[2] +
                        t * (F[3] +
                             t * (F[4] +
                                  t * (F[5] +
                                       t * (F[6] +
                                            t * (F[7] +
                                                 t * (F[8] + t * F[9]))))))))) /
             (F[10] +
              s * (F[11] +
                   s * (F[12] +
                        s * (F[13] +
                             s * (F[14] +
                                  s * (F[15] +
                                       s * (F[16] +
                                            s * (F[17] +
                                                 s * (F[18] + s))))))))) *
             F[19]; // correction to remove bias
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;
        s = one - t;

        fd = (G[0] +
              t * (G[1] +
                   t * (G[2] +
                        t * (G[3] +
                             t * (G[4] +
                                  t * (G[5] +
                                       t * (G[6] +
                                            t * (G[7] +
                                                 t * (G[8] + t * G[9]))))))))) /
             (G[10] +
              s * (G[11] +
                   s * (G[12] +
                        s * (G[13] +
                             s * (G[14] +
                                  s * (G[15] +
                                       s * (G[16] +
                                            s * (G[17] +
                                                 s * (G[18] + s))))))))) *
             G[19]; // correction to remove bias
    }
    else
    {
        w = one / (x * x);
        t = onethousandsixhundred * w;

        fd = x * x * x * x * x * x * x * x * x * x * sqrt(x) * factor *
             (H[0] +
              w * (H[1] +
                   t * (H[2] +
                        t * (H[3] + t * (H[4] + t * (H[5] - t * H[6]))))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 10 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p10(const BS_REAL y)
{
    BS_REAL x, ex, t, s, y2;
    BS_REAL fd = 0.;

    x = -fabs(y);

    constexpr BS_REAL A[6]  = {3628800,
                               1769.11836499240874,
                               2.74114748349296929,
                               0.0153129509121134183,
                               0.000171415624076050339,
                               3.15741345799313491e-6};
    constexpr BS_REAL B[13] = {
        7.23923596166699396e9, 7.46011115766110324e9, 3.51560076972821513e9,
        9.79963332403905796e8, 1.72722667850513422e8, 1.85872779338803730e7,
        980119.027043354796,   5345.1738990257193,    5661.6190861926977,
        2763.5410988949565,    803.48749805350521,    149.54574471403198,
        17.324182954806273};
    constexpr BS_REAL C[6] = {7.25062913034521571e6, 1.20504343904096290e6,
                              59606.1300016688751,   1363.72727447603412,
                              16.4493406684822644,   0.0909090909090909091};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;
        s  = one - t;

        fd = ex *
             (A[0] -
              ex * (A[1] + s * (A[2] + s * (A[3] + s * (A[4] + s * A[5])))));
    }
    else if (x <= zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] + t * (B[4] + t * (B[5] + t * B[6])))))) /
             (B[7] +
              s * (B[8] +
                   s * (B[9] + s * (B[10] + s * (B[11] + s * (B[12] + s))))));
    }
    if (y > zero)
    {
        y2 = y * y;

        fd = fd +
             y * (C[0] +
                  y2 * (C[1] +
                        y2 * (C[2] + y2 * (C[3] + y2 * (C[4] + y2 * C[5])))));
    }

    return fd;
}

/* BS_REAL precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 21/2 */
KOKKOS_INLINE_FUNCTION
BS_REAL FDI_p212(const BS_REAL x)
{
    const BS_REAL factor = 2. / 23.; // = 1/(k+1)
    BS_REAL ex, t, w, s;
    BS_REAL fd = 0.;

    constexpr BS_REAL A[6]  = {1.18994230839622485e7,   4103.25503466860405,
                               5.19741906649655513,     0.0251967194122553158,
                               0.000252954565705636841, 4.24573915094571073e-6};
    constexpr BS_REAL B[13] = {
        2.57000236613575175e10, 2.66775768017824432e10, 1.26823522594150685e10,
        3.57329482988205525e9,  6.38366360660921950e8,  6.99165387034866263e7,
        3.77747237079641640e6,  5829.6137460125218,     6133.1261260721688,
        2967.5333196500928,     853.66205862182849,     156.68183527916256,
        17.793529538010685};
    constexpr BS_REAL C[15] = {
        3.9256463321926022e11, 3.5272008321450459e11, 1.4960338986194452e11,
        4.0840964995951679e10, 8.3525752792093004e9,  1.3652900099806473e9,
        1.6653931710965140e8,  1.13722223824443501e7, 10780.778916243366,
        12211.902429240048,    6810.538241268853,     2444.9651837399508,
        622.4594160475152,     115.21708456233468,    14.653186311439082};
    constexpr BS_REAL D[17] = {
        1.02038438555518024e12, 1.44086845004445038e12, 1.06951376053758491e12,
        5.83349119667798955e11, 2.59184513055987834e11, 9.22344989320630316e10,
        2.46767885435893916e10, 4.50298910194390990e9,  4.40533249962065063e8,
        2627.6481435132529,     3598.485403241250,      2816.3184982451496,
        1591.2052134075520,     693.3250333689612,      234.74251456309597,
        59.835610870372227,     10.529060212037859};
    constexpr BS_REAL E[18] = {
        1.92539440493876721e13, 6.76427905453399646e13, 1.17204818714740525e14,
        1.31312575057612394e14, 1.04191636638567130e14, 5.99784170385635521e13,
        2.45015915838649510e13, 6.53409875052899427e12, 8.86443190814042498e11,
        3208.0283085434934,     3638.3347710577402,     2454.8591095140763,
        1237.5705806010665,     499.5532233797705,      164.10632643205828,
        43.045792336177945,     8.4174363747925276,     0.999999999999999649};
    constexpr BS_REAL F[20] = {
        1.35032839928658828e16, 9.84341732357224541e16, 3.35732954478237964e17,
        7.02136022561204639e17, 9.91445873700748081e17, 9.80394164297009879e17,
        6.79928848238746956e17, 3.19982319837818266e17, 9.33078432378309346e16,
        1.29960135147133318e16, 32616.320007579440,     25996.868562405952,
        13806.569048984017,     5938.718527590630,      2189.8437434853080,
        703.2974320942947,      195.48717008127781,     45.691859712028360,
        8.3721801701244348,     0.999999999999999529};
    constexpr BS_REAL G[20] = {
        1.01458237235662575e19, 9.67883827639700444e19, 4.18744967361252695e20,
        1.07961610867931522e21, 1.83106271883557796e21, 2.12354501642149972e21,
        1.68956127221959607e21, 8.93514699913180709e20, 2.87036545075886579e20,
        4.31795348604572483e19, 32461.378904270757,     24848.502722479329,
        12785.219096902698,     5373.990729362032,      1953.3664643688677,
        624.0842567143388,      174.30734991710066,     41.427358351772922,
        7.8380857354891532,     0.999999999999999445};
    constexpr BS_REAL H[7] = {1,
                              198.625788571923339,
                              11.5426284919561671,
                              0.365993577138676996,
                              0.00572282501335319477,
                              0.0000313848422097767888,
                              1.47431097320479878e-8};
    constexpr BS_REAL half = 0.5;

    if (x < -two)
    {
        ex = exp(x);
        t  = ex * fdi_litconst;
        s  = one - t;

        fd = ex *
             (A[0] -
              ex * (A[1] + s * (A[2] + s * (A[3] + s * (A[4] + s * A[5])))));
    }
    else if (x < zero)
    {
        s = -half * x;
        t = one - s;

        fd = (B[0] +
              t * (B[1] +
                   t * (B[2] +
                        t * (B[3] + t * (B[4] + t * (B[5] + t * B[6])))))) /
             (B[7] +
              s * (B[8] +
                   s * (B[9] + s * (B[10] + s * (B[11] + s * (B[12] + s))))));
    }
    else if (x < two)
    {
        t = half * x;
        s = one - t;

        fd = (C[0] +
              t * (C[1] +
                   t * (C[2] +
                        t * (C[3] +
                             t * (C[4] +
                                  t * (C[5] + t * (C[6] + t * C[7]))))))) /
             (C[8] +
              s * (C[9] +
                   s * (C[10] +
                        s * (C[11] +
                             s * (C[12] + s * (C[13] + s * (C[14] + s)))))));
    }
    else if (x < five)
    {
        t = one_third * (x - two);
        s = one - t;

        fd =
            (D[0] +
             t * (D[1] +
                  t * (D[2] +
                       t * (D[3] +
                            t * (D[4] +
                                 t * (D[5] +
                                      t * (D[6] + t * (D[7] + t * D[8])))))))) /
            (D[9] +
             s * (D[10] +
                  s * (D[11] +
                       s * (D[12] +
                            s * (D[13] +
                                 s * (D[14] +
                                      s * (D[15] + s * (D[16] + s))))))));
    }
    else if (x < ten)
    {
        t = one_fifth * x - one;
        s = one - t;

        fd =
            (E[0] +
             t * (E[1] +
                  t * (E[2] +
                       t * (E[3] +
                            t * (E[4] +
                                 t * (E[5] +
                                      t * (E[6] + t * (E[7] + t * E[8])))))))) /
            (E[9] +
             s * (E[10] +
                  s * (E[11] +
                       s * (E[12] +
                            s * (E[13] +
                                 s * (E[14] +
                                      s * (E[15] + s * (E[16] + s)))))))) *
            E[17]; // correction to remove bias
    }
    else if (x < twenty)
    {
        t = one_tenth * x - one;
        s = one - t;

        fd = (F[0] +
              t * (F[1] +
                   t * (F[2] +
                        t * (F[3] +
                             t * (F[4] +
                                  t * (F[5] +
                                       t * (F[6] +
                                            t * (F[7] +
                                                 t * (F[8] + t * F[9]))))))))) /
             (F[10] +
              s * (F[11] +
                   s * (F[12] +
                        s * (F[13] +
                             s * (F[14] +
                                  s * (F[15] +
                                       s * (F[16] +
                                            s * (F[17] +
                                                 s * (F[18] + s))))))))) *
             F[19]; // correction to remove bias
    }
    else if (x < forty)
    {
        t = one_twentyth * x - one;
        s = one - t;

        fd = (G[0] +
              t * (G[1] +
                   t * (G[2] +
                        t * (G[3] +
                             t * (G[4] +
                                  t * (G[5] +
                                       t * (G[6] +
                                            t * (G[7] +
                                                 t * (G[8] + t * G[9]))))))))) /
             (G[10] +
              s * (G[11] +
                   s * (G[12] +
                        s * (G[13] +
                             s * (G[14] +
                                  s * (G[15] +
                                       s * (G[16] +
                                            s * (G[17] +
                                                 s * (G[18] + s))))))))) *
             G[19]; // correction to remove bias
    }
    else
    {
        w = one / (x * x);
        t = onethousandsixhundred * w;

        fd = x * x * x * x * x * x * x * x * x * x * x * sqrt(x) * factor *
             (H[0] +
              w * (H[1] +
                   t * (H[2] +
                        t * (H[3] + t * (H[4] + t * (H[5] + t * H[6]))))));
    }

    return fd;
}

/*===========================================================================*/

// fermi_distr.c

// Computation of Fermi-Dirac distribution function

/* Inputs:
 * 	e     [MeV] : energy
 * 	temp  [MeV] : temperature
 * 	mu    [MeV] : chemical potential
 */
KOKKOS_INLINE_FUNCTION
BS_REAL FermiDistr(const BS_REAL e, const BS_REAL temp, const BS_REAL mu)
{
    const BS_REAL arg = (e - mu) / temp;
    BS_REAL tmp;

    constexpr BS_REAL one = 1;

    // Handle differently arg>0 and arg<0 cases
    if (arg > zero)
    {
        tmp = SafeExp(-arg);
        return tmp / (tmp + one);
    }
    else
    {
        tmp = SafeExp(arg);
        return one / (tmp + one);
    }
}

/*===========================================================================*/

// Calculates the exponential-dependent factors in the denominator of the NEPS
// kernel
KOKKOS_INLINE_FUNCTION
BS_REAL NEPSExpFunc(BS_REAL x)
{
    const BS_REAL exp_ = SafeExp(-fabs(x));
    const int is_x_neg = signbit(x);

    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;

    return (is_x_neg - (! is_x_neg) * exp_) / (one - exp_ + (x == zero));
}


/*===========================================================================*/

// gamma.c

// Computation of gamma function
// Taken from Numerical Recipes ("numerical.recipes/book/book.html")
KOKKOS_INLINE_FUNCTION
BS_REAL Gammln(const BS_REAL xx)
{
    int j;
    BS_REAL x, tmp, y, ser;

    constexpr BS_REAL cof[14] = {
        57.1562356658629235,     -59.5979603554754912,
        14.1360979747417471,     -0.491913816097620199,
        .339946499848118887e-4,  .465236289270485756e-4,
        -.983744753048795646e-4, .158088703224912494e-3,
        -.210264441724104883e-3, .217439618115212643e-3,
        -.164318106536763890e-3, .844182239838527433e-4,
        -.261908384015814087e-4, .368991826595316234e-5};
    constexpr BS_REAL a = 5.24218750000000000;
    constexpr BS_REAL b = 0.5;
    constexpr BS_REAL c = 0.999999999999997092;
    constexpr BS_REAL d = 2.5066282746310005;

    if (xx <= 0)
    {
        // throw("bad arg in gammln");
        BS_PRINTF("bad arg in gammln");
        // exit(1);
    }

    y = x = xx;
    tmp   = x + a;
    tmp   = (x + b) * log(tmp) - tmp;
    ser   = c;

    for (j = 0; j < 14; j++)
    {
        ser += cof[j] / ++y;
    }

    return tmp + log(d * ser / x);
}

// Compute the Gamma function using Serling's approximation
KOKKOS_INLINE_FUNCTION
BS_REAL GammaStirling(const BS_REAL x)
{
    constexpr BS_REAL zero  = 0;
    constexpr BS_REAL e     = 2.718281828459045235360287471352;
    constexpr BS_REAL twopi = 6.283185307179586;

    BS_ASSERT(x > zero);

    return sqrt(twopi / x) * pow(x / e, x);
}

/*===========================================================================*/

// lambert.c

/*
 * Computation of Lambert function using logarithmic recursive formula
 * in "Loczi, Applied Mathematics and Computation, 433 (2022)"
 */

// Recursive formula for W0 real branch of Lambert function
KOKKOS_INLINE_FUNCTION
BS_REAL W0(const BS_REAL x)
{
    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;
    constexpr BS_REAL e    = 2.718281828459045235360287471352;
    constexpr BS_REAL ie   = one / e;
    BS_REAL beta;

    BS_ASSERT(x > -ie);

    // Initial guess for recursion
    if (x > e)
    {
        BS_REAL log_x = log(x);
        beta          = log_x - log(log_x);
    }
    else if (x > zero)
    {
        beta = x / e;
    }
    else
    { // if(x > -ie)
        BS_REAL tmp_1 = e * x;
        BS_REAL tmp_2 = sqrt(one + tmp_1);
        BS_REAL tmp_3 = one + tmp_2;
        beta          = tmp_1 * log(tmp_3) / (tmp_2 * tmp_3);
    }

    beta = beta / (one + beta) * (one + log(x / beta));
    beta = beta / (one + beta) * (one + log(x / beta));
    beta = beta / (one + beta) * (one + log(x / beta));

    return beta;
}


/*===========================================================================*/

// mnewt.c

// One-dimensional root finding functions

// 1D root-finding parameters
// @TODO: find suitable values for the following parameters
constexpr int ntrial_1d   = 150;     // Maximum allowed number of iterations.
constexpr BS_REAL xacc_1d = 1.0e-07; // Set the accuracy for Newton Raphson

// 1D Newton-Raphson with analytic derivative
KOKKOS_INLINE_FUNCTION
BS_REAL MNewt1d(BS_REAL guess, BS_REAL x1, BS_REAL x2, BS_REAL f0,
                MyFunction* func, MyFunction* dfunc)
{
    // Using a combination of Newton-Raphson and bisection, return the solution
    // of f(x) = f0 bracketed between x1 and x2. The root will be refined
    // until its accuracy is known within xacc. f is a user-supplied struct
    // that returns the function value at the point x. df is a user-supplied
    // struct that returns the value of the function first derivative at the
    // point x.
    BS_REAL xh, xl;
    BS_REAL fl             = func->function(&x1, func->params) - f0;
    BS_REAL fh             = func->function(&x2, func->params) - f0;
    constexpr BS_REAL zero = 0;
    constexpr BS_REAL half = 0.5;
    constexpr BS_REAL two  = 2;

    if ((fl > zero && fh > zero) || (fl < zero && fh < zero))
    {
        BS_PRINTF("xl = %.3e, fl = %.3e\n", x1, fl);
        BS_PRINTF("xh = %.3e, fh = %.3e\n", x2, fh);
        BS_PRINTF("Root must be bracketed in rtsafe");
        exit(1);
        // throw("Root must be bracketed in rtsafe");
    }

    if (fl == zero)
        return x1;
    if (fh == zero)
        return x2;

    if (fl < zero)
    { // Orient the search so that f(xl) < 0.
        xl = x1;
        xh = x2;
    }
    else
    {
        xh = x1;
        xl = x2;
    }

    BS_REAL rts   = guess; // 0.5*(x1+x2);  // Initialize the guess for root,
    BS_REAL dxold = fabs(x2 - x1); // the “stepsize before last,”
    BS_REAL dx    = dxold;         // and the last step.

    BS_REAL f  = func->function(&rts, func->params) - f0;
    BS_REAL df = dfunc->function(&rts, dfunc->params);

    for (int j = 0; j < ntrial_1d; j++)
    { // Loop over allowed iterations.
        if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > zero) ||
            (fabs(two * f) > fabs(dxold * df)))
        { // Bisect if Newton out of range, or not decreasing fast enough.
            dxold = dx;
            dx    = half * (xh - xl);
            rts   = xl + dx;
            if (xl == rts)
                return rts;
        }
        else
        { // Change in root is negligible. Newton step acceptable. Take it.
            dxold        = dx;
            dx           = f / df;
            BS_REAL temp = rts;
            rts -= dx;
            if (temp == rts)
                return rts;
        }

        if (fabs(dx) < xacc_1d)
            return rts; // Convergence criterion.

        f = func->function(&rts, func->params) -
            f0; // The one new function evaluation per iteration.
        df = dfunc->function(&rts, dfunc->params);

        if (f < zero)
        { // Maintain the bracket on the root.
            xl = rts;
        }
        else
        {
            xh = rts;
        }
    }
    // throw("Maximum number of iterations exceeded in rtsafe");
    BS_PRINTF("Maximum number of iterations exceeded in rtsafe");
    exit(1);
}

/*===========================================================================*/

// Two-dimensional root finding functions

// Invert two-dimensional matrix
/*
 * in  : input 2x2 matrix
 * out : output 2x2 matrix
 *
 * Entries:
 * 	0 --> [0,0]
 * 	1 --> [0,1]
 *	2 --> [1,0]
 *	3 --> [1,1]
 */
KOKKOS_INLINE_FUNCTION
void Invert2DMat(BS_REAL* in, BS_REAL* out)
{
    const BS_REAL den = in[0] * in[3] - in[1] * in[2];
    out[0]            = in[3] / den;
    out[1]            = -in[1] / den;
    out[2]            = -in[2] / den;
    out[3]            = in[0] / den;
    return;
}

// Implementation of 2D Newton-Raphson root finding

/*
 * Finding solution of A = C - F[x] = 0 (system of two non-linear coupled
 *equations)
 *
 * Input:
 * 	- x (dim = 2) --> initial guess
 * 	- C (dim = 2) --> inhomogeneous term
 *	- PairF (dim = 2) --> calculated using pointer to fdf(x,C,fvec,fjac)
 *function Output:
 * 	- x (1x2 matrix) --> Newton-Raphson solution after convergence
 *
 * - fvec: A (dim = 2)
 * - fjac: Jacobian matrix of A (dim = 2x2)
 */

// 2D root-finding parameters
// @TODO: do some tests and optimize the following parameters
const int ntrial_2d   = 1000;  // Max number of NR iterations
const BS_REAL tolx_2d = 1.e-5; // Accuracy level on the variable
const BS_REAL tolf_2d = 1.e-7; // Accuracy level on the functions

KOKKOS_INLINE_FUNCTION
void MNewt2d(BS_REAL* x, BS_REAL C[2],
             void (*fdf)(BS_REAL*, BS_REAL*, BS_REAL*, BS_REAL*))
{
    const int n = 2; // Matrix dimension
    int i;
    BS_REAL p[n], fvec[n];
    BS_REAL fjac[n * n]; // Jacobian matrix
    BS_REAL finv[n * n]; // Inverse of Jacobian matrix
    for (int k = 0; k < ntrial_2d; k++)
    {
        fdf(x, C, fvec, fjac);
        BS_REAL errf = 0.0;
        for (i = 0; i < n; i++)
            errf += fabs(fvec[i]);
        if (errf <= tolf_2d)
            return;
        Invert2DMat(fjac, finv);
        for (i = 0; i < n; i++)
            p[i] = -(finv[i * n] * fvec[0] +
                     finv[i * n + 1] * fvec[1]); //-fvec[i];
        BS_REAL errx = 0.0;
        for (i = 0; i < n; i++)
        {
            errx += fabs(p[i]);
            x[i] += p[i];
            // BS_PRINTF("x[%d] = %.3e\n", i, x[i]);
        }
        if (errx <= tolx_2d)
            return;
    }
    return;
}

/*===========================================================================*/

// theta.c

// Implementations of the Heaviside Step function
//
// Continuous approximation is taken from Bierng-Chearl Ahn (2013)

#define a_heaviside 20.

// Implementation with if statement
KOKKOS_INLINE_FUNCTION
BS_REAL HeavisidePiecewise(const BS_REAL x)
{
    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;

    if (x < zero)
    {
        return zero;
    }
    else
    {
        return one;
    }
}

// Continuous approximation with tanh - (Eq.5)
KOKKOS_INLINE_FUNCTION
BS_REAL HeavisideTanhApprox(const BS_REAL x)
{
    constexpr BS_REAL half = 0.5;
    constexpr BS_REAL one  = 1;
    constexpr BS_REAL two  = 2;

    constexpr BS_REAL kBS_aHeaviside = a_heaviside;

    BS_REAL kBS_REAL_MIN;

    if (std::is_same_v<BS_REAL, float>)
    {
        kBS_REAL_MIN = FLT_MIN; // FLT_TRUE_MIN
    }
    else
    {
        kBS_REAL_MIN = DBL_MIN; // DBL_TRUE_MIN
    }

    return half *
           (one + tanh(two * kBS_aHeaviside * x / (fabs(x) + kBS_REAL_MIN)));
}

/*===========================================================================*/

// m1matrix.c
//  \brief Routines for operating on M1Matrix structures
KOKKOS_INLINE_FUNCTION
void InitializeM1MatrixSingleFlavor(M1Matrix* mat, const int n, const int idx)
{
    mat->m1_mat_ab[idx] = (BS_REAL**)malloc(sizeof(BS_REAL*) * 2 * n);
    mat->m1_mat_em[idx] = (BS_REAL**)malloc(sizeof(BS_REAL*) * 2 * n);

    for (int i = 0; i < 2 * n; i++)
    {
        mat->m1_mat_ab[idx][i] = (BS_REAL*)malloc(sizeof(BS_REAL) * 2 * n);
        mat->m1_mat_em[idx][i] = (BS_REAL*)malloc(sizeof(BS_REAL) * 2 * n);

        for (int j = 0; j < 2 * n; j++)
        {
            mat->m1_mat_ab[idx][i][j] = zero;
            mat->m1_mat_em[idx][i][j] = zero;
        }
    }

    return;
}

KOKKOS_INLINE_FUNCTION
void FreeM1MatrixSingleFlavor(M1Matrix* mat, const int n, const int idx)
{
    for (int i = 0; i < 2 * n; i++)
    {
        free(mat->m1_mat_ab[idx][i]);
        free(mat->m1_mat_em[idx][i]);
    }

    free(mat->m1_mat_ab[idx]);
    free(mat->m1_mat_em[idx]);

    return;
}

KOKKOS_INLINE_FUNCTION
void InitializeM1Matrix(M1Matrix* mat, const int n)
{
    for (int idx = 0; idx < total_num_species; idx++)
    {
        InitializeM1MatrixSingleFlavor(mat, n, idx);
    }

    return;
}

KOKKOS_INLINE_FUNCTION
void FreeM1Matrix(M1Matrix* mat, const int n)
{
    for (int idx = 0; idx < total_num_species; idx++)
    {
        FreeM1MatrixSingleFlavor(mat, n, idx);
    }

    return;
}
#endif // BNS_NURATES_SRC_FUNCTIONS_FUNCTIONS_H_
