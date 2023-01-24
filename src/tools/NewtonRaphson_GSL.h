#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multiroots.h>

double x_NRaphson(int max_iter, double guess, gsl_function_fdf FDF) {
	int status;
        int iter = 0;
        double x0, x = guess;
        const gsl_root_fdfsolver_type *T;
        gsl_root_fdfsolver *s;
	T = gsl_root_fdfsolver_newton;
        s = gsl_root_fdfsolver_alloc (T);
        gsl_root_fdfsolver_set (s, &FDF, x);

        do {
		iter++;
		status = gsl_root_fdfsolver_iterate (s);
		x0 = x;
		x = gsl_root_fdfsolver_root (s);
		status = gsl_root_test_delta (x, x0, 0, 1e-5);

		//if (status == GSL_SUCCESS)
			//printf ("Converged:\n");

	}
	while (status == GSL_CONTINUE && iter < max_iter);
	
	return x;
}

void print_state (size_t iter, gsl_multiroot_fdfsolver * s) {
  printf ("iter = %3lu x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
  return;
}


void NRaphson_2D(double x_i[2], double *xf, gsl_multiroot_function_fdf FDF) {
        const gsl_multiroot_fdfsolver_type *T;
        gsl_multiroot_fdfsolver *s;
        gsl_multiroot_function_fdf f = FDF;

        int status;
        size_t i, iter = 0;

        gsl_vector *x = gsl_vector_alloc(f.n);

        gsl_vector_set(x, 0, x_i[0]);
        gsl_vector_set(x, 1, x_i[1]);

        T = gsl_multiroot_fdfsolver_newton;
        s = gsl_multiroot_fdfsolver_alloc(T, f.n);
        gsl_multiroot_fdfsolver_set(s, &f, x);

        //print_state (iter, s);

	do
          {
                iter++;
                status = gsl_multiroot_fdfsolver_iterate (s);

                //print_state (iter, s);

                if (status)
                        break;

                status = gsl_multiroot_test_residual(s->f, 1e-7);
          }
        while (status == GSL_CONTINUE && iter < 1000);

        //printf ("status = %s\n", gsl_strerror (status));

        xf[0] = gsl_vector_get (s->x, 0),
        xf[1] = gsl_vector_get (s->x, 1),
        gsl_multiroot_fdfsolver_free (s);

        gsl_vector_free(x);
        return;
}



void NRaphson_2D_findiff(double x_i[2], double *xf, gsl_multiroot_function F) {
        const gsl_multiroot_fsolver_type *T;
        gsl_multiroot_fsolver *s;
        gsl_multiroot_function f = F;

        int status;
        size_t i, iter = 0;

        gsl_vector *x = gsl_vector_alloc(f.n);

        gsl_vector_set(x, 0, x_i[0]);
        gsl_vector_set(x, 1, x_i[1]);

        T = gsl_multiroot_fsolver_dnewton;
        s = gsl_multiroot_fsolver_alloc(T, f.n);
        gsl_multiroot_fsolver_set(s, &f, x);

        //print_state (iter, s);

        do
          {
                iter++;
                status = gsl_multiroot_fsolver_iterate (s);

                //print_state (iter, s);

                if (status)
                        break;

                status = gsl_multiroot_test_residual(s->f, 1e-7);
          }
        while (status == GSL_CONTINUE && iter < 1000);

        //printf ("status = %s\n", gsl_strerror (status));

        xf[0] = gsl_vector_get (s->x, 0),
        xf[1] = gsl_vector_get (s->x, 1),
        gsl_multiroot_fsolver_free (s);

        gsl_vector_free(x);
        return;
}

