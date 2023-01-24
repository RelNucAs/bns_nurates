#include <stdio.h>
#include <math.h>
#include "tools/integration.hpp"

int main (){
        double a = 0., b = 1.;
	const int ngle = 20.;
	const int ngla = 20.;
        save_gauleg(ngle, a, b);
        save_gaulag(ngla, 0.);
	return 0;
}

