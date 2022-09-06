## Python script to compare emissivities and opacities from Fortran and C++ codes

## Import modules
import numpy as np
import py_parameters as pars
import units
import math
import matplotlib.pyplot as plt

## Define functions for comparison
def code_comparison(tit,arr1,arr2):
    print("\n\n###################################################")
    print("   "+tit+"    ")
    print("###################################################")
    print("Fortran code:")
    print(arr1)
    print("C++ code:")
    print(arr2)
    print("Relative difference:")
    diff = compare_arrays(arr1,arr2)
    print(diff)
    print("Max difference: %.5e" %np.amax(np.nan_to_num(diff)))

def compare_arrays(arr1,arr2):
    arr3 = abs(arr1-arr2)/np.minimum(abs(arr1),abs(arr2))
    arr3[np.where((arr1 == np.nan) | (arr2 == np.nan))] = np.nan
    arr3[np.where((arr1 == 0.) & (arr2 == 0.))] = 0.
    return arr3

def plot_rates(Fdata,Cdata,spec):
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row')

    fig.suptitle(spec['title'])

    for j in range(2):
        axs[0][j].set_xscale('log')

    axs[0][0].set_yscale('log')

    axs[1][0].set_xlabel(r'$r$ [cm]')
    axs[1][1].set_xlabel(r'$\rho$ [g/cm$^3$]')

    axs[0][0].set_ylabel(spec['y1lab'])
    axs[1][0].set_ylabel(spec['y2lab'])

    axs[1][0].set_ylim((-0.0005, 0.0005))

    # Remove (horizontal) space between axes
    fig.subplots_adjust(hspace=0,wspace=0)

    # Plot each graph, and manually set the y tick values
    axs[0][0].plot(r,abs(Fdata),marker='.',label='Fortran')
    axs[0][0].plot(r,abs(Cdata),marker='.',label='C++')

    axs[0][1].plot(d,abs(Fdata),marker='.')
    axs[0][1].plot(d,abs(Cdata),marker='.')

    diff = compare_arrays(Fdata,Cdata)

    axs[1][0].plot(r,diff,marker='.')
    axs[1][0].axhline(y=0., color='k', linestyle='--')

    axs[1][1].plot(d,diff,marker='.')
    axs[1][1].axhline(y=0., color='k', linestyle='--')

    nan_slice = np.argwhere(np.isnan(diff))
    #nan_slice = np.where(diff == np.nan)
    #print(diff)
    #print(nan_slice)
    axs[1][0].vlines(x=r[nan_slice], ymin=axs[1][0].get_ylim()[0], ymax=axs[1][0].get_ylim()[1], linestyles='solid', colors='r')
    axs[1][1].vlines(x=d[nan_slice], ymin=axs[1][1].get_ylim()[0], ymax=axs[1][1].get_ylim()[1], linestyles='solid', colors='r')

    fig.legend(loc='upper right')

    plt.savefig(pars.plotfolder + spec['plotname'])

## Print type of compasion 
if pars.use_dU:
    dU_spec = "corrected by dU    #"
    dU_str  = '_dU'
else:
    dU_spec = "without dU correction #"
    dU_str  = ''
print("#############################################################################################")
print("# Comparing output of Fortran and C++ codes including WM corrections, "+dU_spec)
print("#############################################################################################")
print("")

## Import Fortran data
#f_data = np.nan_to_num(np.loadtxt(pars.f_file,comments='#',unpack=True,max_rows=pars.nrowmax))
f_data = np.loadtxt(pars.f_file,comments='#',unpack=True,max_rows=pars.nrowmax)
r,d,T,ye  = f_data[1:5]
f_rates = f_data[11:24] 

## Import C++ data
#C_data = np.nan_to_num(np.loadtxt(pars.C_file,comments='#',unpack=True,max_rows=nrowmax))
C_data = np.loadtxt(pars.C_file,comments='#',unpack=True,max_rows=pars.nrowmax)
C_rates = C_data[8:21]

# legend =  0: dU,   1: j_nue,  2: 1/l_nue, 3: j_anue,  4: 1/l_anue,
#           5: R (nue), 6: Rbar (anue), 7: B_IS (nu), 8: B_IS (anu),
#           9: R (nu+n), 10: Rbar (anu+n),  11: R (nu+p),   12: Rbar (anu+p)

## Compare dU values
dU   = np.column_stack((f_rates[0],C_rates[0]))
if (np.unique(dU[:,0]-dU[:,1])[0] != 0.):
    print("Different dU in Fortran and C++")
else:
    print("dU values:")
    print(dU[:,0])
print("\n\n")


## Compare WM corrections
R    = np.column_stack((f_rates[5],C_rates[5]))
Rbar = np.column_stack((f_rates[6],C_rates[6]))

Rn    = np.column_stack((f_rates[9], C_rates[9]))
Rbarn = np.column_stack((f_rates[10],C_rates[10]))
Rp    = np.column_stack((f_rates[11],C_rates[11]))
Rbarp = np.column_stack((f_rates[12],C_rates[12]))

if (np.unique(R).size > 2):
    print("R is not unique, R unique size = %d" %np.unique(R).size)
else:
    if ((R[0,0] < 0.) or (R[0,1] < 0.)):
        print("WARNING: negative WM correction")
    else:
        print("R (Fortran) = %.3e, R (C++) = %.3e" %(R[0,0],R[0,1]))

if (np.unique(Rbar).size > 2):
    print("Rbar is not unique, Rbar unique size = %d" %np.unique(Rbar).size)
else:
    if ((R[0,0] < 0.) or (R[0,1] < 0.)):
        print("WARNING: negative WM correction")
    else:
        print("Rbar (Fortran) = %.3e, Rbar (C++) = %.3e" %(Rbar[0,0],Rbar[0,1]))

if (np.unique(Rn).size > 2):
    print("Rn is not unique, Rn unique size = %d" %np.unique(Rn).size)
else:
    if ((Rn[0,0] < 0.) or (Rn[0,1] < 0.)):
        print("WARNING: negative WM correction")
    else:
        print("Rn (Fortran) = %.3e, Rn (C++) = %.3e" %(Rn[0,0],Rn[0,1]))

if (np.unique(Rbarn).size > 2):
    print("Rbarn is not unique, Rbarn unique size = %d" %np.unique(Rbarn).size)
else:
    if ((Rbarn[0,0] < 0.) or (Rbarn[0,1] < 0.)):
        print("WARNING: negative WM correction")
    else:
        print("Rbarn (Fortran) = %.3e, Rbarn (C++) = %.3e" %(Rbarn[0,0],Rbarn[0,1]))

if (np.unique(Rp).size > 2):
    print("Rp is not unique, Rp unique size = %d" %np.unique(Rp).size)
else:
    if ((Rp[0,0] < 0.) or (Rp[0,1] < 0.)):
        print("WARNING: negative WM correction")
    else:
        print("Rp (Fortran) = %.3e, Rp (C++) = %.3e" %(Rp[0,0],Rp[0,1]))

if (np.unique(Rbarp).size > 2):
    print("Rbarp is not unique, Rbarp unique size = %d" %np.unique(Rbarp).size)
else:
    if ((Rbarp[0,0] < 0.) or (Rbarp[0,1] < 0.)):
        print("WARNING: negative WM correction")
    else:
        print("Rbarp (Fortran) = %.3e, Rbarp (C++) = %.3e" %(Rbarp[0,0],Rbarp[0,1]))


## Compare neutrino emissivity and opacity
ab_em = np.column_stack((f_rates[1],f_rates[2],f_rates[3],f_rates[4],C_rates[1],C_rates[2],C_rates[3],C_rates[4]))

tit = 'Electron neutrino emissivity (j_nue)'
idx = 0
code_comparison(tit,ab_em[:,idx],ab_em[:,idx+4]) #Electron neutrino emissivity

tit = 'Electron neutrino opacity (1/l_nue)'
idx = 1
code_comparison(tit,ab_em[:,idx],ab_em[:,idx+4]) #electron neutrino opacity

tit = 'Electron antineutrino emissivity (j_anue)'
idx = 2
code_comparison(tit,ab_em[:,idx],ab_em[:,idx+4]) #electron antineutrino emissivity

tit = 'Electron antineutrino opacity (1/l_anue)'
idx = 3
code_comparison(tit,ab_em[:,idx],ab_em[:,idx+4]) #electron antineutrino opacity


## Compare scattering source term
B_IS  = np.column_stack((f_rates[7],f_rates[8],2*math.pi*C_rates[7],2*math.pi*C_rates[8])) #probably 2pi factor more in Fortran code

tit = 'Neutrino scattering source term (B_IS_nu)'
idx = 0
code_comparison(tit,B_IS[:,idx],B_IS[:,idx+2]) #neutrino scattering

tit = 'Antineutrino scattering source term (B_IS_anu)'
idx = 1
code_comparison(tit,B_IS[:,idx],B_IS[:,idx+2]) #antineutrino scattering


# Make plots for comparison
em_nue = {'title': r'$e^-+p\rightarrow\nu_e+n$:     emissivity',
        'y1lab': r'$j_{\nu_e}$ [1/s]',
        'y2lab': r'$|\Delta j|/j$',
        'plotname': 'em_nue_'+pars.E+dU_str+'.png'}
idx = 0
plot_rates(ab_em[:,idx],ab_em[:,idx+4],em_nue)


ab_nue = {'title': r'$\nu_e+n\rightarrow e^-+p$:  opacity',
        'y1lab': r'$\chi_{\nu_e}$ [1/cm]',
        'y2lab': r'$|\Delta \chi|/\chi$',
        'plotname': 'ab_nue_'+pars.E+dU_str+'.png'}
idx = 1
plot_rates(ab_em[:,idx],ab_em[:,idx+4],ab_nue)


em_anue = {'title': r'$e^++n\rightarrow\bar{\nu}_e+p$:   emissivity',
        'y1lab': r'$j_{\bar{\nu}_e}$ [1/s]',
        'y2lab': r'$|\Delta j|/j$',
        'plotname': 'em_anue_'+pars.E+dU_str+'.png'}
idx = 2
plot_rates(ab_em[:,idx],ab_em[:,idx+4],em_anue)


ab_anue = {'title': r'$\bar{\nu}_e+n\rightarrow e^++n$:   opacity',
        'y1lab': r'$\chi_{\bar{\nu}_e}$ [1/cm]',
        'y2lab': r'$|\Delta \chi|/\chi$',
        'plotname': 'ab_anue_'+pars.E+dU_str+'.png'}
idx = 3
plot_rates(ab_em[:,idx],ab_em[:,idx+4],ab_anue)

B_IS_nu = {'title': r'$\nu+N\rightarrow\nu+N$:  source term',
        'y1lab': r'$B^{\rm IS}_{\nu_e}$ [1/cm]',
        'y2lab':  r'$\Delta|R^{\rm IS}|/R^{\rm IS}$',
        'plotname': 'B_IS_nu_'+pars.E+dU_str+'.png'}
idx = 0
plot_rates(B_IS[:,idx],B_IS[:,idx+2],B_IS_nu)


B_IS_anu = {'title': r'$\bar{\nu}+N\rightarrow\bar{\nu}+N$:  source term',
        'y1lab': r'$B^{\rm IS}_{\bar{\nu}_e}$ [1/cm]',
        'y2lab': r'$\Delta|R^{\rm IS}|/R^{\rm IS}$',
        'plotname': 'B_IS_anu_'+pars.E+dU_str+'.png'}
idx = 1
plot_rates(B_IS[:,idx],B_IS[:,idx+2],B_IS_anu)

