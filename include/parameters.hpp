//Define parameters
#pragma once //compile only once

#include <string>

namespace parameters
{
	//decide to use (1) or not (0) recoil and weak magnetism for neutrino absorption on nucleons 
        const int use_WM_ab = 0;
	
	//decide to use (1) or not (0) recoil and weak magnetism for elastic neutrino scattering on nucleons 
        const int use_WM_sc = 0;

	//decide to use (1) or not (0) deltaU 
        const int use_dU = 1;
	
	//energy bins parameter
	const int ne = 20;
	const double emin = 3.e0;
	const double emax = 3.e2;

	//numerical integration
	const bool save_wghs = true;	
	const double alpha = 0.; //power exponent of Gauss-Laguerre weight
	
	const std::string abs_path = "/home/leonardo/Desktop/PhD_work/BNS_muons/bns_nurates/"; //to be moved somewhere else
}


