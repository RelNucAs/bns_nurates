import numpy as np
import re
import matplotlib.pyplot as plt
from argparse import ArgumentParser 

## Plot emissivity and inverse mean free path along
## 1D CCSN SN profile + comparison with results from
## reference Fortran code

parser = ArgumentParser(description="Inclusion of dU correction")
parser.add_argument('--du', dest='du', default=False, action='store_true')
args = parser.parse_args()

## Inclusion of dU correction
use_dU = args.du
print("use_dU = %r" %use_dU)
print("")

if not use_dU:
    dU_str = ""
else:
    dU_str = "_dU"

print("###############################################");
print("# Comparing spectral rates for beta reactions #")
print("###############################################");

## Function computing the absolute value of the relative difference between two quantities
def rel_diff(y1, y2):
    return abs(y1 - y2) / y2

## Define filename
filename = "../tests/tests_beta/test_beta_1.008e+01" + dU_str + ".txt"

## Get neutrino energy 
omega = np.loadtxt(filename, skiprows=0, max_rows=1, comments=None, dtype=str)[3] # MeV
omega = float(omega)

print("E_nu = %.5e MeV" %omega)
print("")

## Get weak magnetism corrections
w_f  = np.loadtxt(filename, skiprows=18, max_rows=2, comments=None, dtype=str)[:,-1]
w_nue_f  = float(w_f[0])
w_anue_f = float(w_f[1])

w_c  = np.loadtxt(filename, skiprows=21, max_rows=2, comments=None, dtype=str)[:,-1]
w_nue_c  = float(w_f[0])
w_anue_c = float(w_f[1])

## Print weak magnetism corrections
print("W_nue_f  = %lf" %w_nue_f)
print("W_anue_f = %lf" %w_anue_f)
print("")

print("W_nue_c  = %lf" %w_nue_c)
print("W_anue_c = %lf" %w_anue_c)
print("")

## Read profile data
idx = np.loadtxt(filename, comments="#", delimiter=" ", usecols=(0), dtype=int)
r, rho, temp, ye = np.loadtxt(filename, comments="#", delimiter=" ", usecols=(1,2,3,4), unpack=True)
mu_e, mu_hat = np.loadtxt(filename, comments="#", delimiter=" ", usecols=(5,6), unpack=True)
yh, ya, yp, yn = np.loadtxt(filename, comments="#", delimiter=" ", usecols=(7,8,9,10), unpack=True)
du = np.loadtxt(filename, comments="#", delimiter=" ", usecols=(11), unpack=True)
em_nue_f, ab_nue_f, em_anue_f, ab_anue_f = np.loadtxt(filename, comments="#", delimiter=" ", usecols=(12,13,14,15), unpack=True) # Fortran results
em_nue_c, ab_nue_c, em_anue_c, ab_anue_c = np.loadtxt(filename, comments="#", delimiter=" ", usecols=(16,17,18,19), unpack=True) # C results

## Compute relative difference (Fortran vs C)
em_nue_diff  = rel_diff(em_nue_c , em_nue_f)
em_anue_diff = rel_diff(em_anue_c, em_anue_f)

ab_nue_diff  = rel_diff(ab_nue_c , ab_nue_f)
ab_anue_diff = rel_diff(ab_anue_c, ab_anue_f)

## Generate plot
fig, axs = plt.subplots(2, 2, sharex="all", figsize=(8,3), gridspec_kw={'height_ratios': [3, 1]})

plt.suptitle("Comparison of beta rates (dU: %r)" %use_dU)

## Plot of the rates
ax = axs[0][0]
ax.set_title(r"$\nu_e + n \rightarrow e^- + p$ ($E_\nu = %.3e$ MeV)" %omega)
ax.plot(r, em_nue_f, label=r"Fortran - $j_{\nu_{\rm e}}$")
ax.plot(r, ab_nue_f, label=r"Fortran - $\lambda^{-1}_{\nu_{\rm e}}$")
ax.plot(r, em_nue_c, label=r"C - $j_{\nu_{\rm e}}$")
ax.plot(r, ab_nue_c, label=r"C - $\lambda^{-1}_{\nu_{\rm e}}$")
ax.set_ylabel(r"$j~[{\rm s}^{-1}]$  -  $\lambda^{-1}~[{\rm cm}^{-1}]$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()

ax = axs[0][1]
ax.set_title(r"$\bar{\nu}_e + p \rightarrow e^+ + n$ ($E_\nu = %.3e$ MeV)" %omega)
ax.plot(r, em_anue_f, label=r"Fortran - $j_{\bar{\nu}_{\rm e}}$")
ax.plot(r, ab_anue_f, label=r"Fortran - $\lambda^{-1}_{\bar{\nu}_{\rm e}}$")
ax.plot(r, em_anue_c, label=r"C - $j_{\bar{\nu}_{\rm e}}$")
ax.plot(r, ab_anue_c, label=r"C - $\lambda^{-1}_{\bar{\nu}_{\rm e}}$")
ax.set_ylabel(r"$j~[{\rm s}^{-1}]$  -  $\lambda^{-1}~[{\rm cm}^{-1}]$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()

## Plot of the relative difference
ax = axs[1][0]
ax.plot(r, em_nue_diff, label=r"$\Delta j / j$")
ax.plot(r, ab_nue_diff, label=r"$\Delta \lambda^{-1} / \lambda^{-1}$")
ax.set_xlabel(r"$r~[\rm cm]$")
ax.set_yscale("log")
ax.legend()

ax = axs[1][1]
ax.plot(r, em_anue_diff, label=r"$\Delta j / j$")
ax.plot(r, ab_anue_diff, label=r"$\Delta \lambda^{-1} / \lambda^{-1}$")
ax.set_xlabel(r"$r~[\rm cm]$")
ax.set_yscale("log")
ax.legend()

plt.subplots_adjust(hspace=0.0)

plt.show()
