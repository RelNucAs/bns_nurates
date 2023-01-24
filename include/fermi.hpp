// Fermi distribution and related tools
#pragma once //compile only once


//double Fermi(double mu, double E, double T)
//{
	//return 1/(exp((E-mu)/T)+1);
//}


double Fermi(const double mu, const double E, const double T) 
{
  double tmp;
  const double arg = (E-mu) / T;

  if (arg > 0.) {
	  tmp = exp(-arg);
	  return tmp / (tmp+1.);
  } else {
	  tmp = exp(arg);
	  return 1.  / (tmp+1.);
  }
}
