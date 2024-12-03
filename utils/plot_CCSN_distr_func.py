import numpy as np
import re
import matplotlib.pyplot as plt
from argparse import ArgumentParser 

## Plot neutrino equilibrium distribution function
## at a given point of the input 1D CCSN profile

parser = ArgumentParser(description="CCSN profile radius")
parser.add_argument('-r', dest='radius', required=True, type=float)
args = parser.parse_args()

Q = 1.29333236 # nucleon mass difference [MeV]

## Find array index closest to value
def find_index(value, array):
    idx = np.argmin(abs(array-value))
    idx = min(array.size,max(0, idx))
    return idx

## Compute neutrino distribution function at equilibrium
def nu_distribution(omega, temp, mu):
    return 1. / (np.exp((omega - mu) / temp) + 1.)

## Read 1D CCSN input profile
filename = "../inputs/nurates_CCSN/nurates_1.008E+01.txt"

data = np.loadtxt(filename, unpack=True)

r  = data[1] # radius [cm]
d  = data[2] # mass density [g/cm3]
T  = data[3] # temperature [MeV]
ye = data[4] # electron fraction
mu_e   = data[5] # relativistic electron chemical potential [MeV]
mu_hat = data[6] # non-relativistic neutron-proton chemical potential [MeV]

mu_nue  = mu_e - mu_hat - Q
mu_anue = - mu_nue
mu_nux  = np.zeros(mu_nue.size)

## Print radius
idx = find_index(args.radius, r)
print("idx = %d: " %idx + "r = %e cm" %r[idx])

## Thermodynamics conditions
txt = r"$r = %.2e$ cm"            % r[idx] + "\n" \
      r"$\rho = %.2e$ g/cm$^3$"   % d[idx] + "\n" \
      r"$T = %.2e$ MeV"           % T[idx] + "\n" \
      r"$Y_{\rm e} = %.2lf$"      % ye[idx] + "\n" \
      r"$\mu_{\rm e} = %.2e$ MeV" % mu_e[idx] + "\n" \
      r"$\hat{\mu} = %.2e$ MeV"    % mu_hat[idx]

## Neutrino energy array
#omega = np.linspace(0.01, 100., num=1000)
omega = np.logspace(np.log10(0.01), np.log10(1000.), num=1000)

## Generate plot
title = "Neutrino equilibrium distribution function"

fig, ax = plt.subplots(1, 1, sharex="all", num=title)

plt.suptitle(title)

props = dict(boxstyle='round', lw=0.5, facecolor="w", alpha=1.0)
ax.text(0.7, 0.05, txt, fontsize="8.5", transform = ax.transAxes, bbox=props)

ax.plot(omega, nu_distribution(omega, T[idx], mu_nue[idx]) , label=r"$\nu_{\rm e}$ ($\mu = %.2e$ MeV)" %mu_nue[idx])
ax.plot(omega, nu_distribution(omega, T[idx], mu_anue[idx]), label=r"$\bar{\nu}_{\rm e}$ ($\mu = %.2e$ MeV)" %mu_anue[idx])
ax.plot(omega, nu_distribution(omega, T[idx], mu_nux[idx]) , label=r"$\nu_{\rm x} $ ($\mu = %.2e$ MeV)" %mu_nux[idx])

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("Neutrino energy [MeV]")
ax.set_ylabel(r"$f_\nu$")

ax.legend()

plt.show()
