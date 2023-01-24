#pragma once

#include <fstream>

#include "tools/fermi_integrals.h"
#include "tools/nr3.h"
#include "tools/gamma.h"
#include "tools/gauss_wgts.h"
#include "spectral_function.hpp"


double GSL_integration(int n, double a, double b, double alpha, double beta, gsl_function F, const gsl_integration_fixed_type * gauss_type) {
        gsl_integration_fixed_workspace *w;
        double result;
        w = gsl_integration_fixed_alloc(gauss_type, n, a, b, alpha, beta);
        gsl_integration_fixed(&F, &result, w);
        gsl_integration_fixed_free (w);
        return result;
}

