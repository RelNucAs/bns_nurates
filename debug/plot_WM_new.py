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
hnv = -0.5
hna = -0.5*gA
hpv =  0.5-2.*sinsqthetaw
hpa =  0.5*gA

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


def WM_scatt_new(Enu, reacflag):
    [cv,ca,F2] = nucfrmfac(Enu,reacflag) #form factors

    if (reacflag == 1) :
        h0 = hpv*hpv + 3.*hpa*hpa
        h1 = hpv*hpv -    hpa*hpa
    elif (reacflag == 2):
        h0 = hnv*hnv + 3.*hna*hna
        h1 = hnv*hnv -    hna*hna

    ehor = Enu* MeV/(mb*c*c)

    R0 = (cv*cv + 3.*ca*ca + 1.5*ehor*ehor*F2*F2 + 2.*ehor*ehor*cv*F2) / h0
    R1 = (cv*cv -    ca*ca -  2.*ehor*ehor*F2*F2 - 4.*ehor*ehor*cv*F2) / h1

    return [R0,R1]


def scatt_diffCS(Enu, x, reacflag):
    [cv,ca,F2] = nucfrmfac(Enu,reacflag) #form factors

    if (reacflag == 1) :
        h0 = hpv*hpv + 3.*hpa*hpa
        h1 = hpv*hpv -    hpa*hpa
    elif (reacflag == 2):
        h0 = hnv*hnv + 3.*hna*hna
        h1 = hnv*hnv -    hna*hna

    ehor = Enu* MeV/(mb*c*c)

    R0 = (cv*cv + 3.*ca*ca + 1.5*ehor*ehor*F2*F2 + 2.*ehor*ehor*cv*F2) / h0
    R1 = (cv*cv -    ca*ca -  2.*ehor*ehor*F2*F2 - 4.*ehor*ehor*cv*F2) / h1

    return (h0*R0 + h1*R1*x)

## Plot WM as a function of energy

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

[R0_p, R1_p] = WM_scatt_new(E,1)
[R0_n, R1_n] = WM_scatt_new(E,2)

axs[1].plot(E,R0_p,color='tab:red',label=r'$R_0 (\nu p)$')
axs[1].plot(E,R1_p,color='tab:red',ls='--',label=r'$R_1 (\nu p)$')
axs[1].plot(E,R0_n,color='tab:green',label=r'$R_0 (\nu n)$')
axs[1].plot(E,R1_n,color='tab:green',ls='--',label=r'$R_1 (\nu n)$')
#axs[1].axhline(y=0., color='k', linestyle='--')
axs[1].set_ylim((0.4,2.2))
axs[1].legend()

plt.savefig(pars.plotfolder+'WM_vs_energy_new.png',dpi=200)


## Plot differential scattering cross section corrected by WM
fig, ax = plt.subplots(1, sharex='col', sharey='row')

fig.suptitle(r'$\frac{d\sigma}{d\Omega}$ for $\nu+N$ scattering (Horowitz 2002)')

#for j in range(2):
    #axs[j].set_xscale('log')

ax.set_xlabel(r'$E$ [MeV]')
ax.set_ylabel(r'$d\sigma/d\Omega$')

#plot each graph
ax.plot(E,scatt_diffCS(E, 0., 1),label=r'$x=0\,(\nu p)$')
ax.plot(E,scatt_diffCS(E, 0., 2),label=r'$x=0\,(\nu n)$')
ax.plot(E,scatt_diffCS(E, 1., 1),label=r'$x=+1\,(\nu p)$')
ax.plot(E,scatt_diffCS(E, 1., 2),label=r'$x=+1\,(\nu n)$')
ax.plot(E,scatt_diffCS(E, -1., 1),label=r'$x=-1\,(\nu p)$')
ax.plot(E,scatt_diffCS(E, -1., 2),label=r'$x=-1\,(\nu n)$')
ax.legend()

plt.savefig(pars.plotfolder+'diffCS_vs_E.png',dpi=200)


## Plot differential scattering cross section corrected by WM
fig, ax = plt.subplots(1, sharex='col', sharey='row')

fig.suptitle(r'$\frac{d\sigma}{d\Omega}$ for $\nu+N$ scattering (Horowitz 2002)')

#for j in range(2):
    #axs[j].set_xscale('log')

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$d\sigma/d\Omega$')

x = np.linspace(-1,+1,num=1000)

#plot each graph
ax.plot(x,scatt_diffCS(1.,   x, 1),label=r'$E=1\,{\rm MeV}\,(\nu p)$')
ax.plot(x,scatt_diffCS(10.,  x, 1),label=r'$E=10\,{\rm MeV}\,(\nu p)$')
ax.plot(x,scatt_diffCS(100., x, 1),label=r'$E=100\,{\rm MeV}\,(\nu p)$')
ax.plot(x,scatt_diffCS(1.,   x, 2),label=r'$E=1\,{\rm MeV}\,(\nu n)$')
ax.plot(x,scatt_diffCS(10.,  x, 2),label=r'$E=10\,{\rm MeV}\,(\nu n)$')
ax.plot(x,scatt_diffCS(100., x, 2),label=r'$E=100\,{\rm MeV}\,(\nu n)$')
ax.legend()

plt.savefig(pars.plotfolder+'diffCS_vs_x.png',dpi=200)

