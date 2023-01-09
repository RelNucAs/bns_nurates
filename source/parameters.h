//Define parameters
#pragma once //compile only once
#include <string>

struct eos_table_dim {int nne; int nnm; int nt;};

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
	const double eta1_e = -2.3e4;
	const double eta2_e =  2.3e4;
	//muons
	const double eta1_mu = -2.3e4;
        const double eta2_mu = 5.0e4; //1.3e4

	bool HR = true;

	struct eos_table_dim SR_tab = {700 ,  750, 150};
	struct eos_table_dim HR_tab = {1400, 1500, 300};
	
	std::string abs_path = "/home/leonardo/Desktop/PhD_work/BNS_muons/bns_nurates/"; //to be moved somewhere else
}


