#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

// GNU Scientific Library root-finding routines
// Based on GSL 2.7.1

// @TODO: change names of variables and functions following Google C style

/*===========================================================================*/

// One-dimensional root finding functions
// Doc: https://www.gnu.org/software/gsl/doc/html/roots.html

// @TODO: find suitable values for the following parameters (current are taken from GSL examples)
// 1D root-finding parameters
const int max_iter_1d  = 100;
const double epsabs_1d = 0.;
const double epsrel_1d = 1.E-03;

// Print root-finding solver state 
void print_state_1d(size_t iter, const double f, const gsl_root_fdfsolver * s) {
  printf ("iter = %3lu x = % .3f "
          "f(x) = % .3e\n",
          iter,
          s->root,
          f);
  return;
}

// 1D Newton-Raphson with analytic derivative 
double MNewt1d_GSL(const double guess, gsl_function_fdf *fdf) {
  int status, iter = 0;
  double x0, x = guess;
  const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_newton;
  gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc(T);
  gsl_root_fdfsolver_set(s, fdf, x);

  //print_state_1d(iter, fdf->f(s->root,fdf->params), s);
  do {
    iter++;
    status = gsl_root_fdfsolver_iterate(s);
    //print_state_1d(iter, fdf->f(s->root,fdf->params), s);
    x0 = x;
    x = gsl_root_fdfsolver_root(s);
    status = gsl_root_test_delta(x, x0, epsabs_1d, epsrel_1d); // check if |x-x0| < epsabs + epsrel*|x_1|
    //status = gsl_root_test_residual(s->f, epsabs_1d); // check if |f| < epsabs
    //if (status == GSL_SUCCESS) printf ("Converged:\n");
  } while (status == GSL_CONTINUE && iter < max_iter_1d);

  gsl_root_fdfsolver_free(s);  
  return x;
}


/*===========================================================================*/

// Two-dimensional root finding functions
// Doc: https://www.gnu.org/software/gsl/doc/html/multiroots.html

// @TODO: find suitable values for the following parameters (current are taken from GSL examples)
// 2D root-finding parameters
const int max_iter_2d  = 1000;
const double epsabs_2d = 1.E-07;
const double epsrel_2d = 1.E-07;

// Print root-finding solver state 
void print_state_2d(size_t iter, const gsl_multiroot_fdfsolver * s) {
  printf ("iter = %3lu x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
  return;
}

// 2D Newton-Raphson with analytic Jacobian 
int MNewt2d_GSL(double *xi, double *xf, gsl_multiroot_function_fdf *fdf) {
  int status;
  size_t i, iter = 0;
  
  const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_newton;
  gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(T, fdf->n);

  gsl_vector *x = gsl_vector_alloc(fdf->n);
  gsl_vector_set(x, 0, xi[0]);
  gsl_vector_set(x, 1, xi[1]);

  gsl_multiroot_fdfsolver_set(s, fdf, x);

  //print_state_2d(iter, s);

  do {
    iter++;
    status = gsl_multiroot_fdfsolver_iterate(s);
    //print_state_2d(iter, s);
    if (status)  break;
    status = gsl_multiroot_test_residual(s->f, epsabs_2d); // check if \sum |f_i| < epsabs
    //status = gsl_multiroot_test_delta(s->df, s->f, epsabs_2d, epsrel_2d); // check if |dx_i| < epsabs + epsrel*|x_i|
  } while (status == GSL_CONTINUE && iter < max_iter_2d);

  //printf("status = %s\n", gsl_strerror (status));

  xf[0] = gsl_vector_get(s->x, 0),
  xf[1] = gsl_vector_get(s->x, 1),
  gsl_multiroot_fdfsolver_free(s);
  gsl_vector_free(x);
  
  return GSL_SUCCESS;
}


// 2D Newton-Raphson with finite difference Jacobian
int MNewt2d_fd_GSL(double *xi, double *xf, gsl_multiroot_function * f) {
  int status;
  size_t i, iter = 0;
  
  const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_dnewton;
  gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(T, f->n);

  gsl_vector *x = gsl_vector_alloc(f->n);
  gsl_vector_set(x, 0, xi[0]);
  gsl_vector_set(x, 1, xi[1]);

  gsl_multiroot_fsolver_set(s, f, x);

  //print_state_2d(iter, s);

  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate(s);
    //print_state_2d(iter, s);

    if (status) break;
    status = gsl_multiroot_test_residual(s->f, epsabs_2d); // check if \sum |f_i| < epsabs
    //status = gsl_multiroot_test_delta(s->df, s->f, epsabs_2d, epsrel_2d); // check if |dx_i| < epsabs + epsrel*|x_i|

  } while (status == GSL_CONTINUE && iter < max_iter_2d);

  //printf("status = %s\n", gsl_strerror (status));

  xf[0] = gsl_vector_get(s->x, 0),
  xf[1] = gsl_vector_get(s->x, 1),
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free(x);
  
  return GSL_SUCCESS;
}

/*===========================================================================*/
