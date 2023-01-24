// Calculation of single nucleon form factors

#pragma once //compile only once

namespace formfactors
{
		const double lamp = 1.793;
		const double lamn = -1.913;
}


std::tuple<double,double,double> nucfrmfac(const double E, const int reacflag);
