#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

const int nne = 700;
const int nt = 150;

std::string trim(const std::string& str,
                 const std::string& whitespace = " \t")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

void read_electron_table(int *nn1, int *nn2, double *n_arr, double *t_arr, double* eta_arr) {
	string filename = "eos_table/eos_electrons_v2.txt";
	string EOSline;
	int n1, n2;
	int i, j;
	
	ifstream EOStable;
	//EOStable.open(filename,ios::in);  // open
	EOStable.open(filename);  // open
	// open input table
        //std::fstream EOStable("eos_table/eos_electrons_v2.txt"); //read file
        if (!EOStable) {
                cout << "EOS table not found!" << endl; //'\n';
	        std::exit(EXIT_FAILURE);
        }
	
	std::getline(EOStable, EOSline);
	std::stringstream s1(EOSline);
	s1 >> n1;


	std::getline(EOStable, EOSline);
	std::stringstream s2(EOSline);
	s2 >> n2;

        getline(EOStable, EOSline);
	std::stringstream s3(EOSline);
        for (i=0;i<n1;i++) {
        	s3 >> n_arr[i];
        }


        getline(EOStable, EOSline);
	std::stringstream s4(EOSline);
        for (i=0;i<n2;i++) {
        	s4 >> t_arr[i];
        }

	for (i=0;i<n1;i++) {
		getline(EOStable, EOSline);
		std::stringstream s5(EOSline);
		for (j=0;j<n2;j++) {
			s5 >> eta_arr[i*n2+j];
		}
	}
	
        printf("nne = %d, nt = %d\n", n1, n2);

	//std::ifstream EOStable (filename, ios::in | ios::binary);

	// get length of file:
    	//EOStable.seekg (0, EOStable.end);
    	//int length = EOStable.tellg();
    	//EOStable.seekg (0, EOStable.beg);

    	//EOStable.read(reinterpret_cast<char *>(&n1), sizeof(n1));
    	//EOStable.read((char *) &n1, sizeof(n1));
	//EOStable >> n1 >> n2;
	//printf("%d, %d\n", n1, n2);
    	//if (EOStable.tellg() == length) {
	//	std::cout << "all characters read successfully.\n";
    	//} else {
	//	std::cout << "error: not all characters read successfully\n";
	//}
	EOStable.close();
	*nn1 = n1;
	*nn2 = n2;
    	return;
 }

//#include <fstream>
//#include <vector>
//typedef unsigned char BYTE;

//std::vector<BYTE> readFile() {
    // open the file:
    //std::streampos fileSize;
    //string filename = "eos_table/eos_electrons_v2.bin";
    //std::ifstream file(filename, std::ios::binary);

    // get its size:
    //file.seekg(0, std::ios::end);
    //fileSize = file.tellg();
    //file.seekg(0, std::ios::beg);

    // read the data:
    //std::vector<BYTE> fileData(fileSize);
    //file.read((char*) &fileData[0], fileSize);
    //return fileData;
//}

//=======================================================================
//
//     Set the temperature index and interpolation variable
//
//=======================================================================

void sett(double t, double t_arr[], int* idx, double* r) {
	double logt;
	double xmin = log10(t_arr[0]);
	double dx = log10(t_arr[1])-log10(t_arr[0]);
	int it;
        double rt;	

//.....find temperature index.........................................
      	logt = log10(t);
      	it = floor( (logt - xmin)/dx );
      	if (it < 1) {
		it = 0;
	} else if (it >= nt-1) {
        	it = nt - 2;
	}
//.....check the index for non-equidistant grid...........................
      	if (log10(t_arr[it]) > logt) {
        	if (it == 0) {
          		rt = 0.;
		} else {
			it = it - 1;
			rt = (logt - log10(t_arr[it])) / (log10(t_arr[it+1]) - log10(t_arr[it]));
		}
	} else if (log10(t_arr[it+1]) < logt) {
		if (it == nt-2) {
          		rt = 1.;
		} else {
          		it = it + 1;
			rt = (logt - log10(t_arr[it])) / (log10(t_arr[it+1]) - log10(t_arr[it]));
		}
	} else {
        	rt = (logt - log10(t_arr[it])) / (log10(t_arr[it+1]) - log10(t_arr[it]));
	}
      
	if ((rt>1.) || (rt<0.)) {
	       printf("Warning: Something wrong with t index in tinterp!\n");
	       printf("logt = %.3lf, rt = %.3lf\n", logt, rt);
	       printf("it = %d, t[it] = %.3lf, t[it+1] = %.3lf\n", it, t_arr[it], t_arr[it+1]);
               std::exit(EXIT_FAILURE);
      	}

	*idx = it;
	*r   = rt;
	return;
}


//=======================================================================
//
//     Set the density index and interpolation variable
//
//=======================================================================

void setd(double d, double d_arr[], int* idx, double* r) {
	double logd;
	double xmin = log10(d_arr[0]);
	double dx = log10(d_arr[1]) - log10(d_arr[0]);
	int id;
	double rd;

	//.....find density index................................................
      	logd = log10(d);
	id = 1 + floor( (logd - xmin)/dx );
	if (id < 1) {
		id = 0;
	} else if (id >= nne-1) {
		id = nne - 2;
	}
//.....check the index for non-equdistant grid...........................
  	if (log10(d_arr[id]) > logd) {
		if (id == 0) {
        		rd = 0.;
		} else {
			id = id - 1;
			rd = (logd - log10(d_arr[id])) / (log10(d_arr[id+1]) - log10(d_arr[id]));
		}
	} else if (log10(d_arr[id+1]) < logd) {
	        if (id == nne-2) {
          		rd = 1.;
		} else {
			id = id + 1;
	      		rd = (logd - log10(d_arr[id])) / (log10(d_arr[id+1]) - log10(d_arr[id]));
		}
	} else {
        	rd = (logd - log10(d_arr[id])) / (log10(d_arr[id+1]) - log10(d_arr[id]));
	}
      
	if ((rd>1.) || (rd<0.)) {
	       printf("Warning: Something wrong with d index in dinterp!\n");
	       printf("logd = %.3lf, rd = %.3lf\n", logd, rd);
	       printf("id = %d, d[id] = %.3lf, d[id+1] = %.3lf\n", id, d_arr[id], d_arr[id+1]);
               std::exit(EXIT_FAILURE);
      	}
	*idx = id;
	*r   = rd;
	return;
}


//=======================================================================
//
//     EOS: interpolate in eos table with temperature input
//
//=======================================================================


double eos_tintep(double d, double t, double *d_arr, double *t_arr, double *y) {

//     input:
//    d ... density [g/cm^3]
//    t ... temperature [MeV]

        int id, it;
        double rd, rt;
        double yintp;

        // Read table
      //if (teos%initial) then
        //write(6,*) 'Error: eos table unknown in eos_tinterp!'
        //stop
      //endif

        setd(d,d_arr,&id,&rd);
        sett(t,t_arr,&it,&rt);

	printf("%.3lf\n", rt);
//.....interpolate.......................................................
 	 yintp = (1.-rd)*( (1.-rt)*log10(y[nt*id + it])
       	     +                 rt *log10(y[nt*id + it+1]) )
             +       rd *( (1.-rt)*log10(y[nt*(id+1) + it])
       	     +                 rt *log10(y[nt*(id+1) + it+1]) );

        return pow(10.,yintp);
}

