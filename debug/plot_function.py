import numpy as np
import matplotlib.pyplot as plt
import py_parameters as pars

# Profile id 
idx = 48

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

c1 = (GF*GF/pi) / (h*c/(2*pi))**4.

c2 = gV*gV+3*gA*gA

#print(c*c1*c2)

def jnu_function(x, T, mu_e):
    #x = E + Q
    #return c * c1 * c2 * eta_pn * x**2. * (1.-(me/x)**2.)**0.5 / (1. + np.exp((x-mu_e)/T))
    return x**3. * (x+Q)**2. * (1.-(me/(x+Q))**2.)**0.5 / (1. + np.exp((x+Q-mu_e)/T))

def FD_function(x, k, eta):
    return x**k / (np.exp(x-eta)+1.)

def polyFit(x, coeffs):
    dim = coeffs.size
    result = 0.
    for i in range(dim):
        result += coeffs[i] * x**(dim-1-i)
    return result

def compare_arrays(arr1,arr2):
    arr3 = abs(arr1-arr2)/np.minimum(abs(arr1),abs(arr2))
    arr3[np.where((arr1 == np.nan) | (arr2 == np.nan))] = np.nan
    arr3[np.where((arr1 == 0.) & (arr2 == 0.))] = 0.
    return arr3


f_data = np.loadtxt("./nue_em_data.txt") #,unpack=True)
eta_pn = f_data[idx][0]
T      = f_data[idx][1]
mu_e   = f_data[idx][2]

print(f_data[idx][:])

## Plot function as a function of energy
#define energy array
E = np.logspace(np.log10(1.e-3),np.log10(2.e2),num=1000)
E = np.linspace(0.1,200,num=1000)

#create subplot
fig, axs = plt.subplots(2, 1, sharex='col', sharey='row')

fig.suptitle(r'id_r = %d, $\mu_e=%.1lf$, $T=%.2lf$, $\eta_e=%.1lf$' %(idx,mu_e,T,mu_e/T))

axs[0].set_xscale('log')

axs[1].set_xscale('log')
axs[1].set_yscale('log')

axs[1].set_xlabel(r'$E$ [MeV]')
#ax.set_ylabel(r'$j_{\nu}\,(\nu_e)$')

# Remove (horizontal) space between axes
fig.subplots_adjust(hspace=0,wspace=0)

MS = 2

scalingFactor = np.amax(jnu_function(E,T,mu_e)) / np.amax(FD_function((E+Q)/T,5.,mu_e/T))
#plot each graph
for i in range(2):
    axs[i].plot(E,jnu_function(E,T,mu_e),marker='.',ms=MS,label=r'$j_\nu$')
    axs[i].plot(E,FD_function((E+Q)/T,5.,mu_e/T)*scalingFactor,marker='.',ms=MS,label='$F_5(\eta_e)$')
    axs[i].axvline(x=mu_e-Q,color='k',ls='--',lw=0.8,label=r"[$\frac{1}{2}$,1,2]*($\mu_e-Q$)")
    axs[i].axvline(x=0.5*(mu_e-Q),color='k',ls='--',lw=0.8)
    axs[i].axvline(x=2.0*(mu_e-Q),color='k',ls='--',lw=0.8)
plt.legend()

plt.savefig(pars.plotfolder+'jnu_function_idx_'+str(idx)+'.png',dpi=200)


#create subplot
fig, ax = plt.subplots(1, sharex='col', sharey='row')

fig.suptitle(r'Fermi-Dirac integral $F_5(\eta)$')

#ax.set_xscale('log')
#ax.set_yscale('log')

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$F_5(\eta)$')

MS = 2

eta = [0.1,10,100]
t = [5.0359201909e+00,10.,97.5] 
ax.plot(E,FD_function(E,5.,eta[0])/np.amax(FD_function(E,5.,eta[0])),marker='.',ms=MS,color='tab:orange',label='$\eta='+str(eta[0])+'$')
ax.plot(E,FD_function(E,5.,eta[1])/np.amax(FD_function(E,5.,eta[1])),marker='.',ms=MS,color='tab:green' ,label='$\eta='+str(eta[1])+'$')
ax.plot(E,FD_function(E,5.,eta[2])/np.amax(FD_function(E,5.,eta[2])),marker='.',ms=MS,color='tab:blue'  ,label='$\eta='+str(eta[2])+'$')
ax.axvline(x=eta[0],ls='--',lw=0.8,color='tab:orange')
ax.axvline(x=t[0]  ,ls=':' ,lw=0.8,color='tab:orange')
ax.axvline(x=eta[1],ls='--',lw=0.8,color='tab:green')
ax.axvline(x=t[1]  ,ls=':' ,lw=0.8,color='tab:green')
ax.axvline(x=eta[2],ls='--',lw=0.8,color='tab:blue')
ax.axvline(x=t[2]  ,ls=':' ,lw=0.8,color='tab:blue')
ax.legend()

#ax.set_ylim(())
plt.savefig(pars.plotfolder+'FD_function_vs_eta.png',dpi=200)

#create subplot
fig, ax0 = plt.subplots(1, sharex='col', sharey='row')

fig.suptitle(r'Fermi-Dirac integral $F_5(\eta)$')

#ax.set_xscale('log')
#ax.set_yscale('log')

ax0.set_xlabel(r'$x$')
ax0.set_ylabel(r'$F_5(\eta)$')

eta = 100
t = 97.087
ax0.plot(E,FD_function(E,5.,eta)/np.amax(FD_function(E,5.,eta)),marker='.',ms=MS,color='tab:blue',label='$\eta='+str(eta)+'$')
ax0.axvline(x=eta,ls='--',lw=0.8,color='tab:blue')
ax0.axvline(x=t  ,ls=':' ,lw=0.8,color='k')
ax0.legend()

plt.savefig(pars.plotfolder+'FD_function_eta_'+str(eta)+'.png',dpi=200)


data = np.loadtxt("../output/FD_integral_n_32.txt",comments='#',unpack=True)
F_eta = data[0]
F_t   = data[1]

deg = 5
coeffs = np.polyfit(F_eta, F_t, deg)

#create subplot
fig, axs = plt.subplots(3,1, sharex='col', sharey='row')

fig.suptitle(r'Max($F_5(\eta)$)')

#ax1.set_xscale('log')
#axs[0].set_yscale('log')
axs[1].set_yscale('log')
axs[2].set_yscale('log')

axs[1].set_xlabel(r'$\eta$')
#ax1.set_ylabel(r'Max($F_5(\eta)$)')
pos = np.where(F_eta > 0.)


# Remove (horizontal) space between axes
fig.subplots_adjust(hspace=0,wspace=0)

axs[0].axhline(y=5,color='k',ls='--',lw=0.8,label=r"$y=5$")
axs[0].axvline(x=0,color='k',ls='--',lw=0.8)
axs[0].plot(F_eta,F_t,color='tab:orange',label=r'Max($F_5(\eta)$)')
#axs[0].plot(F_eta,polyFit(F_eta, coeffs),marker='.',ms=MS,color='tab:blue',label='Fit') # rcond=None, full=False, w=None, cov=False)
axs[0].plot(F_eta[pos],F_eta[pos],color='tab:green',label=r'$\eta$') # rcond=None, full=False, w=None, cov=False)

axs[0].legend()

x_ref = 15
y_ref = 0.04

#axs[1].axvline(y=y_ref,color='b',ls=':',lw=0.8)
axs[1].axvline(x=x_ref,color='b',ls=':',lw=0.8)
axs[1].axhline(y=5,color='k',ls='--',lw=0.8,label=r"$y=5$")
axs[1].axvline(x=0,color='k',ls='--',lw=0.8)
axs[1].plot(F_eta,F_t,color='tab:orange',label=r'Max($F_5(\eta)$)')
#axs[1].plot(F_eta,polyFit(F_eta, coeffs),marker='.',ms=MS,color='tab:blue',label='Fit') # rcond=None, full=False, w=None, cov=False)
axs[1].plot(F_eta[pos],F_eta[pos],color='tab:green',label=r'$\eta$') # rcond=None, full=False, w=None, cov=False)
axs[1].legend()
axs[1].set_ylim((3.e0,2.e2))

axs[2].axhline(y=y_ref,color='b',ls=':',lw=0.8)
axs[2].axvline(x=x_ref,color='b',ls=':',lw=0.8)
axs[2].axvline(x=0,color='k',ls='--',lw=0.8)
axs[2].plot(F_eta[pos],compare_arrays(F_eta[pos],F_t[pos]),marker='.',ms=MS,color='tab:orange',label=r'Max($F_5(\eta)$)')
axs[2].legend()
axs[2].set_xlabel(r"$\eta$")


xt = axs[2].get_xticks() 
xt = np.append(xt,x_ref)

yt = axs[2].get_yticks() 
yt = np.append(yt,y_ref)

xtl=xt.tolist()
xtl[-1]=str(x_ref)
axs[2].set_xticks(xt)
axs[2].set_xticklabels(xtl)

ytl=yt.tolist()
ytl[-1]=str(y_ref)
axs[2].set_yticks(yt)
axs[2].set_yticklabels(ytl)

axs[0].set_xlim((-50,100))
axs[2].set_ylim((8e-4,2.e0))
plt.savefig(pars.plotfolder+'FD_function_max.png',dpi=200)

