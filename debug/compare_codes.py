## Python script to compare emissivities and opacities from Fortran and C++ codes

import numpy as np
import parameters as pars
import units
import math
import matplotlib.pyplot as plt

# Print 
if pars.use_dU:
    print("Comparing rates corrected by dU")
else:
    print("Comparing rates without dU correction")
#          write(47,77)i,ccsn_r(i),ccsn_d(i),ccsn_T(i),ccsn_ye(i),       &
#     &               mu(je),mu(jhat),y(jh),y(ja),y(jp),y(jn), deltau,   &
#     &               em(ie,len),ab(ie,len)/units%c,                     &
#     &               em(ie,lea),ab(ie,lea)/units%c,                     & 
#     &               R(ie),Rbar(ie),                                    &
#     &               chi(ie,len)/units%c,chi(ie,lea)/units%c,           &
#     &               Rn(ie),Rbarn(ie),Rp(ie),Rbarp(ie)

# Import Fortran data
f_data = np.nan_to_num(np.loadtxt(pars.f_file,comments='#',unpack=True))
r,d,T,ye  = f_data[1:5]
f_rates = f_data[11:24]
        #	fhead = "# id, d [g/cm^3], T [MeV], Ye, mu_e [MeV], mu_hat [MeV], Yp, Yn, dU, em_nue, ab_nue[1/cm], em_anue[1/cm], ab_anue[1/cm], R, Rbar, B_IS_nue[1/cm], B_IS_anue[1/cm], Rn, Rbarn, Rp, Rbarp\n";

# Import C++ data
C_data = np.nan_to_num(np.loadtxt(pars.C_file,comments='#',unpack=True))
C_rates = C_data[8:21]

dU   = np.column_stack((f_rates[0],C_rates[0]))
R    = np.column_stack((f_rates[5],C_rates[5]))
Rbar = np.column_stack((f_rates[6],C_rates[6]))

Rn    = np.column_stack((f_rates[9], C_rates[9]))
Rbarn = np.column_stack((f_rates[10],C_rates[10]))
Rp    = np.column_stack((f_rates[11],C_rates[11]))
Rbarp = np.column_stack((f_rates[12],C_rates[12]))

#print("dU (Fortran) = %.3e, dU (C++) = %.3e" %(dU[0,0],dU[0,1]))
print(dU[:,0]-dU[:,1])

if (np.unique(R).size > 2):
    print("R is not unique, R unique size = %d" %np.unique(R).size)
else:
    print("R (Fortran) = %.3e, R (C++) = %.3e" %(R[0,0],R[0,1]))

if (np.unique(Rbar).size > 2):
    print("Rbar is not unique, Rbar unique size = %d" %np.unique(Rbar).size)
else:
    print("Rbar (Fortran) = %.3e, Rbar (C++) = %.3e" %(Rbar[0,0],Rbar[0,1]))

if (np.unique(Rn).size > 2):
    print("Rn is not unique, Rn unique size = %d" %np.unique(Rn).size)
else:
    print("Rn (Fortran) = %.3e, Rn (C++) = %.3e" %(Rn[0,0],Rn[0,1]))

if (np.unique(Rbarn).size > 2):
    print("Rbarn is not unique, Rbarn unique size = %d" %np.unique(Rbarn).size)
else:
    print("Rbarn (Fortran) = %.3e, Rbarn (C++) = %.3e" %(Rbarn[0,0],Rbarn[0,1]))

if (np.unique(Rp).size > 2):
    print("Rp is not unique, Rp unique size = %d" %np.unique(Rp).size)
else:
    print("Rp (Fortran) = %.3e, Rp (C++) = %.3e" %(Rp[0,0],Rp[0,1]))

if (np.unique(Rbarp).size > 2):
    print("Rbarp is not unique, Rbarp unique size = %d" %np.unique(Rbarp).size)
else:
    print("Rbarp (Fortran) = %.3e, Rbarp (C++) = %.3e" %(Rbarp[0,0],Rbarp[0,1]))

ab_em = np.column_stack((f_rates[1],f_rates[2],f_rates[3],f_rates[4],C_rates[1],C_rates[2],C_rates[3],C_rates[4]))
B_IS  = np.column_stack((f_rates[7],f_rates[8],2*math.pi*C_rates[7],2*math.pi*C_rates[8]))

#for i in range(ab_em.shape[0]):
    #print("%.3e, %.3e" %(ab_em[i,3],ab_em[i,7]))

#for i in range(B_IS.shape[0]):
    #print("%.3e, %.3e" %(B_IS[i,1],B_IS[i,3]))

for i in range(int(ab_em.shape[1]/2)):
    diff = np.nan_to_num(abs((ab_em[:,i]-ab_em[:,i+4])/ab_em[:,i]))
    print(diff)
    print(np.amax(diff))

for i in range(int(B_IS.shape[1]/2)):
    diff = np.nan_to_num(abs((B_IS[:,i]-B_IS[:,i+2])/B_IS[:,i]))
    print(diff)
    print(np.amax(diff[:-1]))

exit()

C_rates[4] = 2*math.pi*C_rates[4] #possibily 2pi factor more in Fortran code 
C_rates[5] = 2*math.pi*C_rates[5]

rel_diff = (f_rates-C_rates)/np.minimum(abs(f_rates),abs(C_rates))
rel_diff = np.where((C_rates == 0.) & (f_rates == 0.), 0., rel_diff)
#print(C_rates[4])
#print(C_rates[5])
rel_diff = np.nan_to_num(rel_diff)
#print([np.amax(rel_diff[i]) for i in range(rel_diff.shape[0])])
exit()

# Make plots for comparison
y1lab = [r'$j_{\nu_e}$ [1/s]', r'$\chi_{\nu_e}$ [1/cm]', r'$j_{\bar{\nu}_e}$ [1/s]', r'$\chi_{\bar{\nu}_e}$ [1/cm]',
         r'$B^{\rm IS}_{\nu_e}$ [1/cm]', r'$B^{\rm IS}_{\bar{\nu}_e}$ [1/cm]']
y2lab = [r'$|\Delta j|/j$', r'$|\Delta \chi|/\chi$', r'$|\Delta j|/j$', r'$|\Delta \chi|/\chi$',
         r'$\Delta|R^{\rm IS}|/R^{\rm IS}$', r'$\Delta|R^{\rm IS}|/R^{\rm IS}$']
tit = [r'$e^-+p\rightarrow\nu_e+n$:     emissivity', r'$\nu_e+n\rightarrow e^-+p$:  opacity', 
        r'$e^++n\rightarrow\bar{\nu}_e+p$:   emissivity', r'$\bar{\nu}_e+n\rightarrow e^++n$:   opacity', 
        r'$\nu+N\rightarrow\nu+N$:  source term', r'$\bar{\nu}+N\rightarrow\bar{\nu}+N$:  source term'] 

for i in range(rel_diff.shape[0]):
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row')

    fig.suptitle(tit[i])

    for j in range(2):
        axs[0][j].set_xscale('log')

    axs[0][0].set_yscale('log')

    axs[1][0].set_xlabel(r'$r$ [cm]')
    axs[1][1].set_xlabel(r'$\rho$ [g/cm$^3$]')

    axs[0][0].set_ylabel(y1lab[i])
    axs[1][0].set_ylabel(y2lab[i])

    axs[1][0].set_ylim((-0.0005, 0.0005))

    # Remove (horizontal) space between axes
    fig.subplots_adjust(hspace=0,wspace=0)

    # Plot each graph, and manually set the y tick values
    axs[0][0].plot(r,abs(f_rates[i]),marker='.',label='Fortran')
    axs[0][0].plot(r,abs(C_rates[i]),marker='.',label='C++')

    axs[0][1].plot(d,abs(f_rates[i]),marker='.')
    axs[0][1].plot(d,abs(C_rates[i]),marker='.')

    axs[1][0].plot(r,rel_diff[i],marker='.')
    axs[1][0].axhline(y=0., color='k', linestyle='--')

    axs[1][1].plot(d,rel_diff[i],marker='.')
    axs[1][1].axhline(y=0., color='k', linestyle='--')

    fig.legend(loc='upper right')

    plt.show()
