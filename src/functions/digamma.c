#include <math.h>
#include <stdio.h>

#include "functions.h"

// TODO: change names of variables and functions following Google C style
 
/*
 * Computation of PairPsi (digamma) function
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

#define EPSILON        2.2204460492503131e-16
#define M_PI_VAL       3.14159265358979323846264338328   /* pi */
#define SQRT_DBL_MIN   1.4916681462400413e-154

typedef struct {
  double val;
  double err;
} SFResult;

struct cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
};

typedef struct cheb_series_struct ChebSeries;

//static
double psics_data[23] = {-.038057080835217922,
                          .491415393029387130, 
                         -.056815747821244730,
                          .008357821225914313,
                         -.001333232857994342,
                          .000220313287069308,
                         -.000037040238178456,
                          .000006283793654854,
                         -.000001071263908506,
                          .000000183128394654,
                         -.000000031353509361,
                          .000000005372808776,
                         -.000000000921168141,
                          .000000000157981265,
                         -.000000000027098646,
                          .000000000004648722,
                         -.000000000000797527,
                          .000000000000136827,
                         -.000000000000023475,
                          .000000000000004027,
                         -.000000000000000691,
                          .000000000000000118,
                         -.000000000000000020};

//static
ChebSeries psi_cs = {psics_data,
                     22,
                     -1, 1,
                     17};

//static
double apsics_data[16] = {-.0204749044678185,
                          -.0101801271534859,
                           .0000559718725387,
                          -.0000012917176570,
                           .0000000572858606,
                          -.0000000038213539,
                           .0000000003397434,
                          -.0000000000374838,
                           .0000000000048990,
                          -.0000000000007344,
                           .0000000000001233,
                          -.0000000000000228,
                           .0000000000000045,
                          -.0000000000000009,
                           .0000000000000002,
                          -.0000000000000000};

//static 
ChebSeries apsi_cs = {apsics_data,
                      15,
                      -1, 1,
                      9};

// Evaluation of the Chebyshev series cs at a given point x
void cheb_eval_e(const ChebSeries * cs,
                 const double x,
                 SFResult * result) {
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  double e = 0.0;

  for (j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
    dd = temp;
  }

  {
    double temp = d;
    d = y*d - dd + 0.5 * cs->c[0];
    e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
  }

  result->val = d;
  result->err = EPSILON * e + fabs(cs->c[cs->order]);

  return; // GSL_SUCCESS;
}


// Evaluation of PairPsi (Digamma) function (result && error)
/* digamma for x both positive and negative; we do both
 * cases here because of the way we use even/odd parts
 * of the function
 */
void psi_x(const double x, SFResult * result) {
  const double y = fabs(x);

  if (x == 0.0 || x == -1.0 || x == -2.0) {
    //result->val = std::nan;
    //result->err = std::nan;
    printf("Error in digamma computation\n");
    exit(EXIT_FAILURE);
  } else if (y >= 2.0) {
    const double t = 8.0/(y*y)-1.0;
    SFResult result_c;
    cheb_eval_e(&apsi_cs, t, &result_c);
    if (x < 0.0) {
      const double s = sin(M_PI_VAL*x);
      const double c = cos(M_PI_VAL*x);
      if(fabs(s) < 2.0*SQRT_DBL_MIN) {
        //result->val = std::nan;
	//result->err = std::nan;
        printf("Error in digamma computation\n");
        exit(EXIT_FAILURE);
      } else {
        result->val  = log(y) - 0.5/x + result_c.val - M_PI_VAL * c/s;
        result->err  = M_PI_VAL*fabs(x)*EPSILON/(s*s);
        result->err += result_c.err;
        result->err += EPSILON * fabs(result->val);
        return; // GSL_SUCCESS;
      }
    } else {
      result->val  = log(y) - 0.5/x + result_c.val;
      result->err  = result_c.err;
      result->err += EPSILON * fabs(result->val);
      return; // GSL_SUCCESS;
    }
  } else { /* -2 < x < 2 */
    SFResult result_c;

    if (x < -1.0) { /* x = -2 + v */
      const double v  = x + 2.0;
      const double t1 = 1.0/x;
      const double t2 = 1.0/(x+1.0);
      const double t3 = 1.0/v;
      cheb_eval_e(&psi_cs, 2.0*v-1.0, &result_c);

      result->val  = -(t1 + t2 + t3) + result_c.val;
      result->err  = EPSILON * (fabs(t1) + fabs(x/(t2*t2)) + fabs(x/(t3*t3)));
      result->err += result_c.err;
      result->err += EPSILON * fabs(result->val);
      return; // GSL_SUCCESS;
    } else if (x < 0.0) { /* x = -1 + v */
      const double v  = x + 1.0;
      const double t1 = 1.0/x;
      const double t2 = 1.0/v;
      cheb_eval_e(&psi_cs, 2.0*v-1.0, &result_c);

      result->val  = -(t1 + t2) + result_c.val;
      result->err  = EPSILON * (fabs(t1) + fabs(x/(t2*t2)));
      result->err += result_c.err;
      result->err += EPSILON * fabs(result->val);
      return; //GSL_SUCCESS;
    } else if (x < 1.0) { /* x = v */
      const double t1 = 1.0/x;
      cheb_eval_e(&psi_cs, 2.0*x-1.0, &result_c);

      result->val  = -t1 + result_c.val;
      result->err  = EPSILON * t1;
      result->err += result_c.err;
      result->err += EPSILON * fabs(result->val);
      return; // GSL_SUCCESS;
    } else { /* x = 1 + v */
      const double v = x - 1.0;
      cheb_eval_e(&psi_cs, 2.0*v-1.0, result);
      return;
    }
  }
}

// Evaluation of PairPsi (Digamma) function (only result)
double SFPsi(const double x) {
  SFResult result;
  psi_x(x, &result);
  return result.val;
}

