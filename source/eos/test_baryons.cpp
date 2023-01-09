#include <array>
#include <vector>
#include "eos_baryons.h"
#include "../parameters.h"
#include "full_eos_module.h" 

int main() {

	std::string EOSpath = parameters::abs_path + "eos_table/dd2_frdm_eos_shen98format_v1.02.tab";

	//eos_type dd2_eos;
	//dd2_eos = eos_readin();
	std::array<double,26> eos_ext;
	double rho = 1.e10;
	double temp = 1.e1;
	
	eos_ext = eos_extended(rho, temp, 0.4, 0.3, 0.2);

	for (int i=0;i<26;i++) printf("eos_extended[%d] = %.10e\n", i, eos_ext[i]);

	std::array<double,53> eos_compl;
	eos_compl = eos_full(rho, temp, 0.4, 0.3);

	for (int i=0;i<53;i++) printf("eos_full[%d] = %.10e\n", i, eos_compl[i]);


        //for (int i=0;i<t_entries;i++) {
        //        for (int j=0;j<y_entries;j++) {
        //            for(int k=0;k<nb_entries;k++) {
        //                for (int l=0;l<entries;l++) {
        //                        printf("%.10e\n", dd2_eos[i][j][k][l]);
        //                }
        //            }
        //        }
        //}

	return 0;
}
