#include <math.h>

#include "functions.h"

// Cut-off parameter for exponential evaluation
// Warning: this value was tuned on the previous version of the code
//          (with C++ math functions) 
const double fExpUppLim = 700.;          // possible alternative choice: const*DBL_MAX_EXP
const double fExpLowLim = -fExpUppLim;   //                              const*DBL_MIN_EXP

// Safe exp function to avoid underflow/overflow
double SafeExp (const double x) {
	return exp(fmin(fmax(x,fExpLowLim),fExpUppLim));
}






