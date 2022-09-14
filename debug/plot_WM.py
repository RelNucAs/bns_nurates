import numpy as np
import matplotlib.pyplot as plt
import py_parameters as pars

## Define constants
MeV = 1.602176634e-6 #conversion MeV to CGS
c = 2.997924562e+10 #cm/s
mb = 1.674e-24 #g
gA = 1.23
gV = 1.
gS = 0.
sinsqthetaw = 0.2325

## Function for computing nuclear form factors
def nucfrmfac(E, reacflag):
    lamp = 1.793
    lamn = -1.913

    ehor = E * MeV/(mb*c*c)

    tau = 0.5*ehor*ehor/(1.+ehor)
    eta = 1./(1.+5.6*tau)
    G = 1./pow(1.+4.97*tau,2)
    Fp1 = (1.+tau*(1.+lamp))*G/(1.+tau)
    Fp2 = lamp*G/(1.+tau)
    Fn1 = tau*lamn*(1.-eta)*G/(1.+tau)
    Fn2 = lamn*(1.+tau*eta)*G/(1.+tau)

    #form factors
    if (reacflag == 1):
        cv  = (0.5-2.*sinsqthetaw)*Fp1 - 0.5*Fn1
        ca  = 0.5*(gA-gS)/pow(1.+3.53*tau,2)
        F2  = (0.5-2.*sinsqthetaw)*Fp2 - 0.5*Fn2
    elif (reacflag == 2):
        cv = (0.5-2.*sinsqthetaw)*Fn1 - 0.5*Fp1
        ca = -0.5*(gA+gS)/pow(1.+3.53*tau,2)
        F2 = (0.5-2.*sinsqthetaw)*Fn2 - 0.5*Fp2
    elif (reacflag == 3):
        cv = Fp1 - Fn1
        ca = gA/pow(1.+3.53*tau,2)
        F2 = Fp2 - Fn2
    else:
        print("Error: reacflag out of range")
        return
    
    return [cv,ca,F2]


## Function for computing WM correction for nue emission/absorption
def WM_nue_abs(e_nu):
    [cv,ca,F2] = nucfrmfac(e_nu,3) #form factors

    ehor = e_nu * MeV/(mb*c*c)
    tmp1 = cv*cv*(1.+4.*ehor+16./3.*ehor*ehor) + 3.*ca*ca*pow(1.+4./3.*ehor,2.) + 8./3.*cv*F2*ehor*ehor + 5./3.*ehor*ehor*(1.+2./5.*ehor)*F2*F2
    tmp2 = 4.*(cv+F2)*ca*ehor*(1.+4./3.*ehor)
    tmp3 = (cv*cv+3.0*ca*ca)*pow(1.+2.*ehor,3.)
    
    return (tmp1+tmp2)/tmp3
        


## Function for computing WM correction for anue emission/absorption
def WM_anue_abs(e_nu):
    [cv,ca,F2] = nucfrmfac(e_nu,3) #form factors

    ehor = e_nu * MeV/(mb*c*c)
    tmp1 = cv*cv*(1.+4.*ehor+16./3.*ehor*ehor) + 3.*ca*ca*pow(1.+4./3.*ehor,2.) + 8./3.*cv*F2*ehor*ehor + 5./3.*ehor*ehor*(1.+2./5.*ehor)*F2*F2
    tmp2 = 4.*(cv+F2)*ca*ehor*(1.+4./3.*ehor)
    tmp3 = (cv*cv+3.0*ca*ca)*pow(1.+2.*ehor,3.)

    return (tmp1-tmp2)/tmp3
    


## Function for computing WM correction for (a)nu scattering on nucleons
def WM_scatt(Enu, reacflag):
    [cv,ca,F2] = nucfrmfac(Enu,reacflag) #form factors
    x = 0. #assume an average angle x=0

    ehor = Enu* MeV/(mb*c*c)
    tmp = (4.*ca*(cv+F2)) / (cv*cv*(1.+x) + ca*ca*(3.-x))
    R    = (1.+(tmp-3.) *ehor*(1.-x))
    Rbar = (1.+(-tmp-3.)*ehor*(1.-x))
                
    return  [R,Rbar]


## Plot WM as a function of energyi

#define energy array
E = np.logspace(np.log10(1.e-3),np.log10(300.),num=1000)

#create subplot
fig, axs = plt.subplots(2, sharex='col', sharey='row')

fig.suptitle('WM corrections')

#for j in range(2):
    #axs[j].set_xscale('log')

axs[1].set_xlabel(r'$E$ [MeV]')

axs[0].set_ylabel('Emission/absoprtion')
axs[1].set_ylabel('Scattering')

#remove (horizontal) space between axes
fig.subplots_adjust(hspace=0,wspace=0)

#plot each graph
axs[0].plot(E,WM_nue_abs(E),label=r'$R$') #,abs(Fdata),marker='.',label='Fortran')
axs[0].plot(E,WM_anue_abs(E),label=r'$\bar{R}$') #,abs(Fdata),marker='.',label='Fortran')
axs[0].legend()

[Rp, Rbarp] = WM_scatt(E,1)
[Rn, Rbarn] = WM_scatt(E,2)

axs[1].plot(E,Rp,color='tab:red',label=r'$R (\nu p)$')
axs[1].plot(E,Rbarp,color='tab:red',ls='--',label=r'$\bar{R} (\bar{\nu} p)$')
axs[1].plot(E,Rn,color='tab:green',label=r'$R (\nu n)$')
axs[1].plot(E,Rbarn,color='tab:green',ls='--',label=r'$\bar{R} (\bar{\nu} n)$')
axs[1].axhline(y=0., color='k', linestyle='--')
axs[1].set_ylim((-0.1,1.1))
axs[1].legend()

plt.savefig(pars.plotfolder+'WM_vs_energy.png',dpi=200)



