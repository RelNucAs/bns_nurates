//Unit conversions

namespace constants
{
	const double eV = 1.6*10e-12

	//Global Constants (CGS)

	const double h_bar = 1.05457266*(10e-27)
	const double c = 29979245800
	const double Gf = 1.436*(10e-49)
	const double me = 9.10938356*(10e-28)
	const double mn = 1.6749286*(10e-24)
	const double mp = 1.6726231*(10e-24)
	const double pi = 3.14159265358979323846
	const double kb = 1.380658*(10e-16)

	//Simplying constants
	const double e_rm = me*c*c
	const double n_rm = mn*c*c
	const double p_rm = mp*c*c
	const double sigma_o = (4*(Gf*Gf)*(e_rm*e_rm))/(pi*pow(h_bar*c,4))
	printf("Sigma_0",sigma_o);
	const double gA = 1.23
	const delta_np = 1.29332*eV*10e6
}