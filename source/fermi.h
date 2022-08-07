// Fermi distribution and related tools

#pragma once //compile only once

namespace fermi
{

  //double Fermi(double mu, double E, double T)
	//{
  		//return 1/(exp((E-mu)/T)+1);
	//}


  double Fermi(double mu, double E, double T) 
  {
	  double tmp;
	  double fermi;
	  double arg = (E-mu) / T;

	  if (arg > 0.) {
		  tmp = exp(-arg);
		  fermi = tmp / (tmp+1.);
	  } else {
		  tmp = exp(arg);
		  fermi = 1.  / (tmp+1.);
	  }

	  return fermi;
  }

}
