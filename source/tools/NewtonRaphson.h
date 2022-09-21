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
