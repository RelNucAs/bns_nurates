//Define parameters
#pragma once //compile only once
#include <string>

struct eos_table_dim {int nne; int nnm; int nt;};
struct eos_table_limits {const double den_min; const double den_max; const double y_min; const double y_max; const double t_min; const double t_max; const double eta_min; const double eta_max;};

namespace parameters
{
	//decide to use (1) or not (0) recoil and weak magnetism for neutrino absorption on nucleons 
        const int use_WM_ab = 0;
	
	//decide to use (1) or not (0) recoil and weak magnetism for elastic neutrino scattering on nucleons 
        const int use_WM_sc = 0;

	//decide to use (1) or not (0) deltaU 
        const int use_dU = 0;
	
	//energy bins parameter
	const int ne = 20;
	const double emin = 3.e0;
	const double emax = 3.e2;

	//numerical integration
	bool save_wghs = true;	
	const double alpha = 0.; //power exponent of Gauss-Laguerre weight
	const int ngle = 20;
	const int ngla = 20;

	//fermionic eos
	//electrons
	struct eos_table_limits e_tab_lim  = {5.e+03, 2.e+16, 5.e-03, 5.5e-01, 5.e-02, 1.05e+02, -2.3e+04, 2.3e+04}; 
	//muons
	struct eos_table_limits mu_tab_lim = {5.e-08, 5.e-01, 1.e+08,  2.e+16, 5.e-02, 1.05e+02, -2.3e+04, 5.0e+04}; //1.3e+04
					   //den_min, den_max, y_min,   y_max,  t_min,    t_max,  eta_min,  eta_max


	bool HR = false;

	struct eos_table_dim SR_tab = {700 ,  750, 150};
	struct eos_table_dim HR_tab = {1400, 1500, 300};
	
	std::string abs_path = "/home/leonardo/Desktop/PhD_work/BNS_muons/bns_nurates/"; //to be moved somewhere else
}


