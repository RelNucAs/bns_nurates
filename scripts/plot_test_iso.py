import numpy as np
import re
import matplotlib.pyplot as plt

## Plot total isoenergetic scattering kernel along
## 1D CCSN SN profile + comparison with results from
## reference Fortran code

print("############################################");
print("# Comparing isoenergetic scattering kernel #")
print("############################################");

## Function computing the absolute value of the relative difference between two quantities
def rel_diff(y1, y2):
    return abs(y1 - y2) / y2

## Define filename
filename = "../tests/tests_iso/test_iso_1.008e+01.txt"

## Get neutrino energy 
omega = np.loadtxt(filename, skiprows=0, max_rows=1, comments=None, dtype=str)[3] # MeV
omega = float(omega)

print("E_nu = %.5e MeV" %omega)
print("")

## Get weak magnetism corrections
w_f  = np.loadtxt(filename, skiprows=17, max_rows=4, comments=None, dtype=str)[:,-1]
w0_n_f = float(w_f[0])
w1_n_f = float(w_f[1])
w0_p_f = float(w_f[2])
w1_p_f = float(w_f[3])

w_c  = np.loadtxt(filename, skiprows=22, max_rows=4, comments=None, dtype=str)[:,-1]
w0_n_c = float(w_c[0])
w1_n_c = float(w_c[1])
w0_p_c = float(w_c[2])
w1_p_c = float(w_c[3])

## Print weak magnetism corrections
print("W0_n_f = %lf" %w0_n_f)
print("W1_n_f = %lf" %w1_n_f)
print("W0_p_f = %lf" %w0_p_f)
print("W1_p_f = %lf" %w1_p_f)
print("")

print("W0_n_c = %lf" %w0_n_c)
print("W1_n_c = %lf" %w1_n_c)
print("W0_p_c = %lf" %w0_p_c)
print("W1_p_c = %lf" %w1_p_c)
print("")

## Read profile data
idx = np.loadtxt(filename, comments="#", delimiter=" ", usecols=(0), dtype=int)
r, rho, temp, ye = np.loadtxt(filename, comments="#", delimiter=" ", usecols=(1,2,3,4), unpack=True)
mu_e, mu_hat = np.loadtxt(filename, comments="#", delimiter=" ", usecols=(5,6), unpack=True)
yh, ya, yp, yn = np.loadtxt(filename, comments="#", delimiter=" ", usecols=(7,8,9,10), unpack=True)
r_iso_f = np.loadtxt(filename, comments="#", delimiter=" ", usecols=(11), unpack=True) # Fortran results
r_iso_c = np.loadtxt(filename, comments="#", delimiter=" ", usecols=(12), unpack=True) # C results

## Compute relative difference (Fortran vs C)
r_iso_diff  = rel_diff(r_iso_c , r_iso_f)

## Generate plot
fig, axs = plt.subplots(2, 1, sharex="all", figsize=(6,6), gridspec_kw={'height_ratios': [3, 1]})

plt.suptitle(r"Comparison of isoenergetic scattering kernel")

## Plot of the rates
ax = axs[0]
ax.set_title(r"$E_\nu = %.3e$ MeV" %omega)
ax.plot(r, r_iso_f, label=r"Fortran")
ax.plot(r, r_iso_c, label=r"C")
ax.set_ylabel(r"$R_{\rm iso}~[{\rm cm}^{-1}]$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()

## Plot of the relative difference
ax = axs[1]
ax.plot(r, r_iso_diff)
ax.set_xlabel(r"$r~[\rm cm]$")
ax.set_ylabel(r"$|\Delta R_{\rm iso}| / R_{\rm iso}$")
ax.set_yscale("log")

plt.subplots_adjust(hspace=0.0)

plt.show()