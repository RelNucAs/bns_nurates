#include <math.h>
#include <assert.h>
#include "functions.h"

// Computation of gamma function
// Taken from Numerical Recipes ("numerical.recipes/book/book.html")
double Gammln(const double xx) {
  int j;
  double x, tmp, y, ser;
  static const double cof[14]={57.1562356658629235,
                              -59.5979603554754912,
			       14.1360979747417471,
			       -0.491913816097620199,
			         .339946499848118887e-4,
			         .465236289270485756e-4,
			        -.983744753048795646e-4,
			         .158088703224912494e-3,
			        -.210264441724104883e-3,
			         .217439618115212643e-3,
			        -.164318106536763890e-3,
			         .844182239838527433e-4,
			        -.261908384015814087e-4,
			         .368991826595316234e-5};
  if (xx <= 0) throw("bad arg in gammln");
  y=x=xx;
  tmp = x+5.24218750000000000;
  tmp = (x+0.5)*log(tmp)-tmp;
  ser = 0.999999999999997092;
  for (j=0;j<14;j++) ser += cof[j]/++y;
  return tmp+log(2.5066282746310005*ser/x);
}

// Compute the Gamma function using Serling's approximation
double GammaStirling(const double x)
{
    static const double e  = 2.718281828459045235360287471352;
    static const double twopi  = 6.283185307179586;

    assert(x > 0.);

    return sqrt(twopi / x) * pow(x / e, x);
}