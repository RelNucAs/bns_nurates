void read_eos_table(double ne[], double T[], double eta[]) {
	int i = 0;
	double filename = "..write filename here..";
	std::fstream fin(filename);

	if (!fin) {
		cout << "EOS table not found" << "\n";
	}

	while (!fin.fail()){
		// read data from input table
                //std::stringstream ss(finLine);

                ss >> ne[i] >> T[i] >> eta[i];
                //printf("ne = %.3e, T = %.3e, eta = %.3e\n", ne[i], T[i], eta[i]);
	}
	return;
}

//=======================================================================
//
//     EOS: interpolate in eos table with temperature input
//
//=======================================================================

double eos_tintep(double d, double T) {

//     input:
//    d ... density [g/cm^3]
//    t ... temperature [MeV]
	
	double y[3][3];
	int id, it;
      	double rd, rt;
      	double y_intp;

	// Read table
      //if (teos%initial) then
        //write(6,*) 'Error: eos table unknown in eos_tinterp!'
        //stop
      //endif

      	setd(d,&id,&rd);
	sett(t,&it,&rt);

//.....interpolate.......................................................
  yintp =  (1.-rd)*( (1.-rt)*y(id ,it)
       +                 rt *y(id ,it+1) )
       +       rd *( (1.-rt)*y(id+1 ,it)
       +                 rt *y(id+1,it+1) );
    	
 	return y_intp;
}


//=======================================================================
//
//     Set the temperature index and interpolation variable
//
//=======================================================================

void sett(double t, double t_arr[], int* it, double* rt) {
	double logt;
	double xmin = t_arr[0];
	double dx = t_arr[1]-t_arr[0];

//.....find temperature index.........................................
      	logt = log10(t);
      	it = 1 + floor( (logt - xmin/dx) )
      	if (it < 1) {
		it = 1;
	} else if (it >= nt) {
        	it = nt - 1;
	}
//.....check the index for non-equidistant grid...........................
      	if (t_arr[it] > logt) {
        	if (it == 1) {
          		rt = 0.
		} else {
			it = it - 1;
			rt = (logt - t_arr[it]) / (t_arr[it+1] - t_arr[it]);
		}
	} else if (t_arr[it+1] < logt) {
		if (it == nt-1) {
          		rt = 1.;
		} else {
          		it = it + 1;
			rt = (logt - t_arr[it]) / (t_arr[it+1] - t_arr[it]);
		}
	} else {
        	rt = (logt - t_arr[it]) / (t_arr[it+1] - t_arr[it]);
	}
      
	if ((rt>1.) || (rt<0.)) {
	       printf("Warning: Something wrong with t index in tinterp!\n");
	       printf("logt = %.3lf, rt = %.3lf\n", logt, rt);
	       printf("it = %d, t[it] = %.3lf, t[it+1] = %.3lf\n", it, t_arr[it], t_arr[it+1]);
               std::exit(EXIT_FAILURE);
      	}
	return;
}


//=======================================================================
//
//     Set the density index and interpolation variable
//
//=======================================================================

void setd(double d, double d_arr[], int* id, double* rd) {
	double logd;
	double xmin = d_arr[0];
	double dx = d_arr[1] - d_arr[0];

	//.....find density index................................................
      	logd = log10(d);
	id = 1 + floor( (logd - xmin/dx) );
	if (id < 1) {
		id = 1;
	} else if (id >= nne) {
		id = nne - 1;
	}
//.....check the index for non-equdistant grid...........................
  	if (d_arr[id] > logd) {
		if (id == 1) {
        		rd = 0.;
		} else {
			id = id - 1;
			rd = (logd - d_arr[id]) / (d_arr[id+1] - d_arr[id]);
		}
	} else if (d_arr[id+1] < logd) {
	        if (id == nne-1) {
          		rd = 1.;
		} else {
			id = id + 1;
	      		rd = (logd - d_arr[id]) / (d_arr[id+1] - d_arr[id]);
		}
	} else {
        	rd = (logd - d_arr[id]) / (d_arr[id+1] - d_arr[id]);
	}
      
	if ((rd>1.) || (rd<0.)) {
	       printf("Warning: Something wrong with d index in dinterp!\n");
	       printf("logd = %.3lf, rd = %.3lf\n", logd, rd);
	       printf("id = %d, d[id] = %.3lf, d[id+1] = %.3lf\n", id, t_arr[id], t_arr[id+1]);
               std::exit(EXIT_FAILURE);
      	}
	return;
}

