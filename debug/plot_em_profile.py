import numpy as np
import matplotlib.pyplot as plt
import py_parameters as pars

## Define constants
Q = 1.2935 #MeV
GF = 8.957e-44 #MeV*cm^3
pi = 3.1415926535898
h = 4.1356943e-21 #Mev*s
MeV = 1.602176634e-6 #conversion MeV to CGS
c = 2.997924562e+10 #cm/s
me = 0.510998928 #MeV
mb = 1.674e-24 #g
gA = 1.23
gV = 1.

def compare_arrays(arr1,arr2):
    arr3 = abs(arr1-arr2)/np.minimum(abs(arr1),abs(arr2))
    arr3[np.where((arr1 == np.nan) | (arr2 == np.nan))] = np.nan
    arr3[np.where((arr1 == 0.) & (arr2 == 0.))] = 0.
    return arr3

def plot_emissivity(idx): #idx: profile index
    data_0 = np.loadtxt("../output/j_function.txt",       skiprows=idx, max_rows= 1)
    data_1 = np.loadtxt("../output/j_function_WM.txt",    skiprows=idx, max_rows= 1)
    data_2 = np.loadtxt("../output/j_function_WM_dU.txt", skiprows=idx, max_rows= 1)
   

    therm = np.loadtxt("../input/nurates_1.008E+01.txt",comments='#',unpack=True)
    mu_e = therm[5][idx]
    T    = therm[3][idx]

    ## Plot function as a function of energy
    #define energy array
    x0 = np.log10(0.1)
    x1 = np.log10(300.)
    nslice = 1001
    dx = (x1-x0) / nslice
    x = np.array([10**(x0+i*dx) for i in range(nslice)])
    #E = np.linspace(0.1,200,num=1000)

    #create subplot
    fig, axs = plt.subplots(1, sharex='col', sharey='row')

    fig.suptitle(r'id_r = %d, $\mu_e=%.1lf$, $T=%.2lf$, $\eta_e=%.1lf$' %(idx,mu_e,T,mu_e/T))

    #axs.set_xscale('log')
    #axs.set_yscale('log')

    axs.set_xlabel(r'$E$ [MeV]')
    axs.set_ylabel(r'$j_{\nu}$')

    MS = 2

    max_0 = np.amax(data_0)
    max_1 = np.amax(data_1)
    max_2 = np.amax(data_2)
   
    x_0 = x[np.argmax(data_0)]
    x_1 = x[np.argmax(data_1)]
    x_2 = x[np.argmax(data_2)]
    
    #plot each graph
    axs.plot(x,data_0/max_0,label=r'$j_\nu$',color='tab:blue')
    axs.plot(x,data_1/max_1,label=r'$j_\nu$',color='tab:green')
    axs.plot(x,data_2/max_2,label=r'$j_\nu$',color='tab:orange')
    axs.axvline(x=x_0,ls='--',lw=0.8,color='tab:blue')
    axs.axvline(x=x_1,ls='--',lw=0.8,color='tab:green')
    axs.axvline(x=x_2,ls='--',lw=0.8,color='tab:orange')
    axs.axvline(x=max(T*5.,mu_e-Q),ls='--',lw=0.8,color='k')
    axs.legend()

    print(x_0)
    axs.set_xlim((x_0*0.5,x_0*3.0))

    plt.savefig(pars.plotfolder+'em_profile/idx_'+str(idx)+'.png',dpi=100)
    plt.close()
    return


for idx in range(102):
    plot_emissivity(idx)
