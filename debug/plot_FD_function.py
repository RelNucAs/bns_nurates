import numpy as np
import matplotlib.pyplot as plt
import py_parameters as pars
import argparse

parser = argparse.ArgumentParser(description="FD integrand parameters")
parser.add_argument("-eta", dest="parsEta", required=True, type=float)
parser.add_argument("-k"  , dest="parsK" , required=True, type=float)
args = parser.parse_args() #pass arguments


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

k = args.parsK
eta = args.parsEta

## Plot function as a function of x
#define x array
#x = np.logspace(np.log10(1.e-3),np.log10(2.e2),num=1000)
x = np.linspace(0.1,200,num=10000)

#create subplot
fig, ax = plt.subplots(1, sharex='col', sharey='row')

fig.suptitle(r'NR Fermi-Dirac integrand $F_%d(\eta=%.1lf)$' %(k,eta))

#ax.set_xscale('log')
#ax.set_yscale('log')

ax.set_xlabel(r'x')
ax.set_ylabel(r'$F_%d(%.1lf)$' %(k,eta))

MS = 2

FD_max = np.amax(FD_function(x,k,eta))
x_max = x[np.argmax(FD_function(x,k,eta))]

ax.plot(x,FD_function(x,k,eta)/FD_max,marker='.',ms=MS) #,label='$F_5(\eta_e)$')
ax.axvline(x=x_max,color='tab:blue',ls='--',lw=0.8, label=r'$F^{\rm max}_5(\eta_e)=%.2lf$' %x_max)
if (eta > 0.):
    ax.axvline(x=eta,color='k',ls='--',lw=0.8,label=r"$\eta=%.1lf$" %eta)
plt.legend()

if (eta > 10.):
    ax.set_xlim((0.2*eta,2.0*x_max))
else:
    ax.set_xlim((0.,25.))

plt.show()


