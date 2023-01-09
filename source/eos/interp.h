#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../parameters.h"

using namespace parameters;
using namespace std;

struct EOSeta {
        //double d[nne]; //density
        //double t[nt]; //temperature
        //double eta[nne*nt]; //relativistic parameter
	int nn; //number of density points
	int nt; //number of temperature points
	std::vector<double> d; //density
	std::vector<double> t; //temperature
	std::vector<double> eta; //relativistic parameter
};

struct EOScomplete {
	int nn; //number of density points
	int nt; //number of temperature points
	std::vector<double> d; //density
	std::vector<double> t; //temperature
	std::vector<double> mu_part; //chemical potential of particles
	std::vector<double> n_part; //number density of particles 
	std::vector<double> n_antipart; //number density of antiparticles 
	std::vector<double> P_part; //pressure of particles 
	std::vector<double> P_antipart; //pressure of antiparticles 
	std::vector<double> e_part; //internal energy of particles 
	std::vector<double> e_antipart; //internal energy of antiparticles
	std::vector<double> s_part; //entropy of particles
	std::vector<double> s_antipart; //entropy of antiparticles
};

void reset_string(std::stringstream ss) {
	ss.str("");
	ss.clear(); // Clear state flags.
	return;
}

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

struct EOSeta read_eta_table(const int sp) {
	std::string filename; //"eos_table/eos_electrons_v2.txt";
	std::string res;
	if (parameters::HR == true) {
		res = "_HR";
	} else {
		res = "";
	}

	if (sp == 1) {
		filename = abs_path + "eos_table/electrons/eos_electrons_leo"+res+".txt"; //_HR.txt";
	} else if (sp == 2) {
		filename = abs_path + "eos_table/muons/eos_muons_leo"+res+".txt";
	} else {
		printf("Number of lepton species: ERROR\n");
		std::exit(EXIT_FAILURE);
	}
		
	std::string EOSline;
	struct EOSeta output;
	int n1, n2;
	int i, j;
	double number;

	std::ifstream EOStable;
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
	output.nn = n1;
	//if (n1 != nne) {
	//	cout << "Number of density entries in the table is not what expected!" << endl;
	//        std::exit(EXIT_FAILURE);
	//}

	std::getline(EOStable, EOSline);
	std::stringstream s2(EOSline);
	s2 >> n2;
	output.nt = n2;
	//if (n2 != nt) {
	//	cout << "Number of temperature entries in the table is not what expected!" << endl;
	//        std::exit(EXIT_FAILURE);
	//}

        getline(EOStable, EOSline);
	std::stringstream s3(EOSline);
        for (i=0;i<n1;i++) {
        	s3 >> number;
		output.d.push_back(number);
        }


        getline(EOStable, EOSline);
	std::stringstream s4(EOSline);
        for (i=0;i<n2;i++) {
        	s4 >> number;
		output.t.push_back(number);
        	//s4 >> output.t[i];
        }

	for (i=0;i<n1;i++) {
		getline(EOStable, EOSline);
		std::stringstream s5(EOSline);
		for (j=0;j<n2;j++) {
			s5 >> number;
			output.eta.push_back(number);
			//s5 >> output.eta[i*n2+j];
		}
	}

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
    	return output;
 }

struct EOScomplete read_total_table(const int sp) {
	std::string filename; //"eos_table/eos_electrons_complete.txt";
        std::string res;
        if (parameters::HR == true) {
                res = "_HR";
        } else {
                res = "";
        }

	if (sp == 1) {
		filename = abs_path + "eos_table/electrons/eos_electrons_complete_leo"+res+".txt";
	} else if (sp == 2) {
		filename = abs_path + "eos_table/muons/eos_muons_complete_leo"+res+".txt";
	} else {
		printf("Number of lepton species: ERROR\n");
		std::exit(EXIT_FAILURE);
	}
		
	string EOSline;
	struct EOScomplete output;
	int n1, n2;
	int i, j;
	double number;

	ifstream EOStable;
	EOStable.open(filename);  // open
        //std::fstream EOStable("eos_table/eos_electrons_v2.txt"); //read file
        if (!EOStable) {
                cout << "EOS table not found!" << endl; //'\n';
	        std::exit(EXIT_FAILURE);
        }

        std::getline(EOStable, EOSline);
        std::stringstream s1(EOSline);
        s1 >> n1;
	output.nn = n1;
        //if (n1 != nne) {
        //        cout << "Number of density entries in the table is not what expected!" << endl;
        //        std::exit(EXIT_FAILURE);
        //}

        std::getline(EOStable, EOSline);
        std::stringstream s2(EOSline);
        s2 >> n2;
	output.nt = n2;
        //if (n2 != nt) {
        //        cout << "Number of temperature entries in the table is not what expected!" << endl;
        //        std::exit(EXIT_FAILURE);
        //}

        //for (j=0;j<4;j++) {
	//	getline(EOStable, EOSline);
        //	std::stringstream s3(EOSline);
        //	for (i=0;i<150;i++) {
        //        	s3 >> output.d[j*150+i];
        //	}
	//}

	//reset_string(s3);
	getline(EOStable, EOSline);
	std::stringstream s3(EOSline);
        for (i=0;i<n1;i++) { //for (i=600;i<700;i++)
		s3 >> number;
		output.d.push_back(number);
        	//s3 >> output.d[i];
        }
	

        getline(EOStable, EOSline);
	std::stringstream s4(EOSline);
        for (i=0;i<n2;i++) {
		s4 >> number;
		output.t.push_back(number);
        	//s4 >> output.t[i];
        }

	for (i=0;i<n1;i++) {
		getline(EOStable, EOSline);
		std::stringstream s5(EOSline);
		for (j=0;j<n2;j++) {
			s5 >> number;
			output.mu_part.push_back(number);
			//s5 >> output.mu_part[i*n2+j];
		}
	}
	
        for (i=0;i<n1;i++) {
                getline(EOStable, EOSline);
                std::stringstream s5(EOSline);
                for (j=0;j<n2;j++) {
			s5 >> number;
			output.P_part.push_back(number);
                        //s5 >> output.P_part[i*n2+j];
                }
        }

        for (i=0;i<n1;i++) {
                getline(EOStable, EOSline);
                std::stringstream s5(EOSline);
                for (j=0;j<n2;j++) {
			s5 >> number;
			output.P_antipart.push_back(number);
                        //s5 >> output.P_antipart[i*n2+j];
                }
        }

        for (i=0;i<n1;i++) {
                getline(EOStable, EOSline);
                std::stringstream s5(EOSline);
                for (j=0;j<n2;j++) {
			s5 >> number;
			output.e_part.push_back(number);
                        //s5 >> output.e_part[i*n2+j];
                }
        }

        for (i=0;i<n1;i++) {
                getline(EOStable, EOSline);
                std::stringstream s5(EOSline);
                for (j=0;j<n2;j++) {
			s5 >> number;
			output.e_antipart.push_back(number);
                        //s5 >> output.e_antipart[i*n2+j];
                }
        }

        for (i=0;i<n1;i++) {
                getline(EOStable, EOSline);
                std::stringstream s5(EOSline);
                for (j=0;j<n2;j++) {
			s5 >> number;
			output.s_part.push_back(number);
                        //s5 >> output.s_part[i*n2+j];
                }
        }

        for (i=0;i<n1;i++) {
                getline(EOStable, EOSline);
                std::stringstream s5(EOSline);
                for (j=0;j<n2;j++) {
			s5 >> number;
			output.s_antipart.push_back(number);
                        //s5 >> output.s_antipart[i*n2+j];
                }
        }

	EOStable.close();
    	return output;
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

void sett(double t, std::vector<double> t_arr, int* idx, double* r) {
	double logt, t1, t2;
	double xmin = t_arr[0];
	double dx = log10(t_arr[1])-log10(t_arr[0]);
	int it, itmp;
        double rt;	

	int nt = t_arr.size();

//.....find temperature index.........................................
      	logt = log10(t);
      	it = floor( (logt - log10(xmin))/dx );
      	if (it < 1) {
		it = 0;
	} else if (it >= nt-1) {
        	it = nt - 2;
	}

	itmp = it;

	t1 = log10(t_arr[it]);
	t2 = log10(t_arr[it+1]);
	
//.....check the index for non-equidistant grid...........................
      	if (t1 > logt) {
        	if (it == 0) {
          		rt = 0.;
		} else {
			it = it - 1;
			t1 = log10(t_arr[it]);
			t2 = log10(t_arr[it+1]);
			rt = (logt - t1) / (t2 - t1);
		}
	} else if (t2 < logt) {
		if (it == nt-2) {
          		rt = 1.;
		} else {
          		it = it + 1;
			t1 = log10(t_arr[it]);
			t2 = log10(t_arr[it+1]);
			rt = (logt - t1) / (t2 - t1);
		}
	} else if (t1 == logt) {
		rt = 0.;
	} else if (t2 == logt) {
		rt = 1.;		
	} else {
        	rt = (logt - t1) / (t2 - t1);
	}
      
	if ((rt>1.) || (rt<0.)) {
	       printf("Warning: Something wrong with t index in tinterp!\n");
	       printf("logt = %.3lf, rt = %.3lf\n", logt, rt);
	       printf("it = %d, itmp = %d, t[it] = %.3lf, t[it+1] = %.3lf\n", it, itmp, t_arr[it], t_arr[it+1]);
               std::exit(EXIT_FAILURE);
      	}

	*idx = it;
	*r   = rt;
	//printf("logt = %.3lf, rt = %.3lf\n", logt, rt);
	//printf("it = %d, t[it] = %.3lf, t[it+1] = %.3lf\n", it, t_arr[it], t_arr[it+1]);
	//printf("idx = %d, r = %.3lf\n", it, rt);
	return;
}


//=======================================================================
//
//     Set the density index and interpolation variable
//
//=======================================================================

void setd(double d, std::vector<double> d_arr, int* idx, double* r) {
	double logd, d1, d2;
	double xmin = d_arr[0];
	double dx = log10(d_arr[1]) - log10(d_arr[0]);
	int id, itmp;
	double rd;

	int nne = d_arr.size();

	//.....find density index................................................
      	logd = log10(d); //d
	id = floor( (logd - log10(xmin))/dx );
	if (id < 1) {
		id = 0;
	} else if (id >= nne-1) {
		id = nne - 2;
	}

	itmp = id;

	d1 = log10(d_arr[id]);
	d2 = log10(d_arr[id+1]);

//.....check the index for non-equdistant grid...........................
  	if (d1 > logd) {
		if (id == 0) {
        		rd = 0.;
		} else {
			id = id - 1;
			if (log10(d_arr[id]) > logd) {
				id = id - 1;
			}
			d1 = log10(d_arr[id]);
			d2 = log10(d_arr[id+1]);
			rd = (logd - d1) / (d2 - d1);
		}
	} else if (d2 < logd) {
	        if (id == nne-2) {
          		rd = 1.;
		} else {
			id = id + 1;
			if (log10(d_arr[id+1]) < logd) {
				id = id + 1;
			}
			d1 = log10(d_arr[id]);
			d2 = log10(d_arr[id+1]);
	      		rd = (logd - d1) / (d2 - d1);
		}
	} else {
        	rd = (logd - d1) / (d2 - d1);
	}
      
	if ((rd>1.) || (rd<0.)) {
	       printf("Warning: Something wrong with d index in dinterp!\n");
	       printf("logd = %.15e, rd = %.3lf\n", logd, rd);
	       printf("logd = %.15e, rd = %.3lf\n", d, rd);
	       printf("id = %d, itmp = %d, d[id] = %.15e, d[id+1] = %.15e\n", id, itmp, d_arr[id], d_arr[id+1]);
	       printf("id = %d, itmp = %d, d[id] = %.15lf, d[id+1] = %.15lf\n", id, itmp, d1, d2);
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


double eos_tintep(double d, double t, std::vector<double> d_arr, std::vector<double> t_arr, std::vector<double> y) {

//     input:
//    d ... density [g/cm^3]
//    t ... temperature [MeV]

        int id, it;
        double rd, rt;
        double yintp;

        int nne = d_arr.size();
        int nt  = t_arr.size();

        // Read table
      //if (teos%initial) then
        //write(6,*) 'Error: eos table unknown in eos_tinterp!'
        //stop
      //endif

        setd(d,d_arr,&id,&rd);
        sett(t,t_arr,&it,&rt);


//.....interpolate.......................................................
 	 yintp = (1.-rd)*( (1.-rt)*y[nt*id + it]
       	     +                 rt *y[nt*id + it+1] )
             +       rd *( (1.-rt)*y[nt*(id+1) + it]
       	     +                 rt *y[nt*(id+1) + it+1] );

        return yintp;
}

