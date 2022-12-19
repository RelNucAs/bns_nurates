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
	std::vector<double> d; //density
	std::vector<double> t; //temperature
	std::vector<double> eta; //relativistic parameter
};

struct EOScomplete {
	std::vector<double> d; //density
	std::vector<double> t; //temperature
	std::vector<double> mu_part; //chemical potential of particles
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

struct EOSeta read_eta_table() {
	//std::string filename = "eos_table/eos_electrons_v2.txt";
	//if (HR == true) {
		std::string filename = "../eos_table/eos_electrons_v2.txt"; //leo_HR.txt";
	//} else {
	//	std::string filename = "../eos_table/eos_electrons_leo.txt";
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
	if (n1 != nne) {
		cout << "Number of density entries in the table is not what expected!" << endl;
	        std::exit(EXIT_FAILURE);
	}

	std::getline(EOStable, EOSline);
	std::stringstream s2(EOSline);
	s2 >> n2;
	if (n2 != nt) {
		cout << "Number of temperature entries in the table is not what expected!" << endl;
	        std::exit(EXIT_FAILURE);
	}

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

struct EOScomplete read_total_table() {
	//string filename = "eos_table/eos_electrons_complete.txt";
	string filename = "../eos_table/eos_electrons_complete_leo_HR.txt";
	//if (HR == true) {
	//	filename = "../eos_table/eos_electrons_complete_leo_HR.txt";
	//} else {
	//	filename = "../eos_table/eos_electrons_complete_leo.txt";
	//}

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
        if (n1 != nne) {
                cout << "Number of density entries in the table is not what expected!" << endl;
                std::exit(EXIT_FAILURE);
        }

        std::getline(EOStable, EOSline);
        std::stringstream s2(EOSline);
        s2 >> n2;
        if (n2 != nt) {
                cout << "Number of temperature entries in the table is not what expected!" << endl;
                std::exit(EXIT_FAILURE);
        }

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
	double logt;
	double xmin = t_arr[0];
	double dx = log10(t_arr[1])-log10(t_arr[0]);
	int it;
        double rt;	

//.....find temperature index.........................................
      	logt = t;
      	it = floor( (log10(logt) - log10(xmin))/dx );
      	if (it < 1) {
		it = 0;
	} else if (it >= nt-1) {
        	it = nt - 2;
	}
//.....check the index for non-equidistant grid...........................
      	if (t_arr[it] > logt) {
        	if (it == 0) {
          		rt = 0.;
		} else {
			it = it - 1;
			rt = (logt - t_arr[it]) / (t_arr[it+1] - t_arr[it]);
		}
	} else if (t_arr[it+1] < logt) {
		if (it == nt-2) {
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
	double logd;
	double xmin = d_arr[0];
	double dx = log10(d_arr[1]) - log10(d_arr[0]);
	int id;
	double rd;

	//.....find density index................................................
      	logd = d; //log10(d);
	id = 1 + floor( (log10(logd) - log10(xmin))/dx );
	if (id < 1) {
		id = 0;
	} else if (id >= nne-1) {
		id = nne - 2;
	}
//.....check the index for non-equdistant grid...........................
  	if (d_arr[id] > logd) {
		if (id == 0) {
        		rd = 0.;
		} else {
			id = id - 1;
			rd = (logd - d_arr[id]) / (d_arr[id+1] - d_arr[id]);
		}
	} else if (d_arr[id+1] < logd) {
	        if (id == nne-2) {
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


double eos_tintep(double d, double t, std::vector<double> d_arr, std::vector<double> t_arr, std::vector<double> y) {

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


//.....interpolate.......................................................
 	 yintp = (1.-rd)*( (1.-rt)*y[nt*id + it]
       	     +                 rt *y[nt*id + it+1] )
             +       rd *( (1.-rt)*y[nt*(id+1) + it]
       	     +                 rt *y[nt*(id+1) + it+1] );

        return yintp;
}

