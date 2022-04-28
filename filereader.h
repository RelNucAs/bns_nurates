#pragma once //compile only once

#include <iostream> //Input&Output stream model based library
#include <fstream> // both read and write from/to files

namespace filereader
{
	ifstream input_file("./Input/nurates_1.008E+01.txt");
	if(input_file.is_open())
	{
		int zone;
	  	double r, rho, temp, ye, mu_e, mu_np, yh, ya, yp, yn;

	  	while (input_file >> zone >> r >> rho >> temp >> ye >> mu_e >> mu_np >> yh >> ya >> yp >> yn)
		{
      		cout << rho << " " << temp << " " << ye << " " << mu_e << " " << " " << mu_np << " " << yp << " " << yn << '\n' << endl;
     	}
     input_file.close();
    }
}