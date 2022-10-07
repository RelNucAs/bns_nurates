//Define parameters
#pragma once //compile only once

namespace parameters
{
	//decide to use (1) or not (0) recoil and weak magnetism for neutrino absorption on nucleons 
        const int use_WM_ab = 1;
	
	//decide to use (1) or not (0) recoil and weak magnetism for elastic neutrino scattering on nucleons 
        const int use_WM_sc = 1;

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
	const double eta1 = -5.e3;
	const double eta2 =  5.e3;
}
