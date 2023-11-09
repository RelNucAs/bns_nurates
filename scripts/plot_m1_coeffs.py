import numpy as np
import re
import matplotlib.pyplot as plt
from argparse import ArgumentParser 

## Compare M1 source coefficients with output of
## WeakRates library along 1D CCSN input profile

## Use standard or modified output data
parser = ArgumentParser(description="Comparison with modified data")
parser.add_argument('--mod', dest='mod', default=False, action='store_true')
args = parser.parse_args()

if args.mod:
    print("Using modified data to align with WeakRates definitions")
    print("")
    
id_std = 0
id_nob = 1

def read_WeakRates_data(filename):
    foldername = "../inputs/wr_standard/"
    data_1 = np.loadtxt(foldername + filename, comments="#", delimiter=" ", unpack=True)

    foldername = "../inputs/wr_noblock/"
    data_2 = np.loadtxt(foldername + filename, comments="#", delimiter=" ", unpack=True)

    return [data_1, data_2]

def read_nurates_data(filename):
    foldername = "../tests/tests_opacities_m1/output/"
    if args.mod: filename = filename.replace(".txt","_mod.txt")
    data = np.loadtxt(foldername + filename, comments="#", delimiter=" ", unpack=True)
    return data

## Read WeakRates data
eta_pair_WR = read_WeakRates_data("eta_pair.dat")
eta_brem_WR = read_WeakRates_data("eta_Bremsstrahlung.dat")
eta_plas_WR = read_WeakRates_data("eta_plasmon.dat")
eta_beta_WR = read_WeakRates_data("eta_betadecay.dat")


k_abs_WR = read_WeakRates_data("k_absorption_np.dat")
k_sct_WR = read_WeakRates_data("k_scattering_np.dat")

r = eta_beta_WR[id_std][0]


## Read bns_nurates data
beta_nurates = read_nurates_data("m1_opacities_abs_em.txt")
iso_nurates  = read_nurates_data("m1_opacities_isoscatt.txt")
pair_nurates = read_nurates_data("m1_opacities_pair.txt")
brem_nurates = read_nurates_data("m1_opacities_brem.txt")


################################
## Plot emission coefficients ##
################################

## Beta decay

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

axs[0][0].set_xscale("log")
plt.subplots_adjust(hspace=0.)

# Suptitle
plt.suptitle("Emissivity coefficient for beta reactions")

# Title
axs[0][0].set_title(r"$\nu_{\rm e}$")
axs[0][1].set_title(r"$\bar{\nu}_{\rm e}$")
axs[0][2].set_title(r"$\nu_{\rm x}$")

# Ylabels
axs[0][0].set_ylabel(r"$\eta~[{\rm MeV}\,{\rm cm}^{-3}\,{\rm s}^{-1}]$")
axs[1][0].set_ylabel("Ratio")

ratio_std = beta_nurates[5:7] / eta_beta_WR[id_std][4:6]
ratio_nob = beta_nurates[5:7] / eta_beta_WR[id_nob][4:6]

for idx in range(2):
    # Eta coefficient
    axs[0][idx].plot(r, eta_beta_WR[id_std][idx+4], color = "tab:blue", ls = "-", label="WR - std")
    axs[0][idx].plot(r, eta_beta_WR[id_nob][idx+4], color = "tab:blue", ls = "--", label="WR - no bl")

    axs[0][idx].plot(r, beta_nurates[idx + 5], color = "tab:orange", label="NR")

    # Ratio
    axs[1][idx].plot(r, ratio_std[idx], color = "tab:blue", ls = "-" , label="WR - std")
    axs[1][idx].plot(r, ratio_nob[idx], color = "tab:blue", ls = "--", label="WR - no bl")

for ax in axs[0]:
    ax.set_yscale("log")
    ax.legend()

for ax in axs[:][1]:
    ax.axhline(1., color = "k", ls = "--", zorder=1)
    ax.set_xlabel(r"$r~[{\rm cm}]$")
    ax.set_yscale("log")

##########################################################

## Pair process

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

axs[0][0].set_xscale("log")
plt.subplots_adjust(hspace=0.)

# Suptitle
plt.suptitle(r"Emissivity coefficient for $e^\pm$ annihilation")

# Title
axs[0][0].set_title(r"$\nu_{\rm e}$")
axs[0][1].set_title(r"$\bar{\nu}_{\rm e}$")
axs[0][2].set_title(r"$\nu_{\rm x}$")

# Ylabels
axs[0][0].set_ylabel(r"$\eta~[{\rm MeV}\,{\rm cm}^{-3}\,{\rm s}^{-1}]$")
axs[1][0].set_ylabel("Ratio")

ratio_std = pair_nurates[5:8] / eta_pair_WR[id_std][4:7]
ratio_nob = pair_nurates[5:8] / eta_pair_WR[id_nob][4:7]

for idx in range(3):
    # Eta coefficient
    axs[0][idx].plot(r, eta_pair_WR[id_std][idx+4], color = "tab:blue", ls = "-", label="WR - std")
    axs[0][idx].plot(r, eta_pair_WR[id_nob][idx+4], color = "tab:blue", ls = "--", label="WR - no bl")

    axs[0][idx].plot(r, pair_nurates[idx + 5], color = "tab:orange", label="NR")

    # Ratio
    axs[1][idx].plot(r, ratio_std[idx], color = "tab:blue", ls = "-" , label="WR - std")
    axs[1][idx].plot(r, ratio_nob[idx], color = "tab:blue", ls = "--", label="WR - no bl")

for ax in axs[0]:
    ax.set_yscale("log")
    ax.legend()

for ax in axs[:][1]:
    ax.axhline(1., color = "k", ls = "--", zorder=1)
    ax.set_xlabel(r"$r~[{\rm cm}]$")
    ax.set_yscale("log")

##########################################################

## Bremsstrahlung process

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

axs[0][0].set_xscale("log")
plt.subplots_adjust(hspace=0.)

# Suptitle
plt.suptitle("Emissivity coefficient for NN bremsstrahlung")

# Title
axs[0][0].set_title(r"$\nu_{\rm e}$")
axs[0][1].set_title(r"$\bar{\nu}_{\rm e}$")
axs[0][2].set_title(r"$\nu_{\rm x}$")

# Ylabels
axs[0][0].set_ylabel(r"$\eta~[{\rm MeV}\,{\rm cm}^{-3}\,{\rm s}^{-1}]$")
axs[1][0].set_ylabel("Ratio")

ratio_std = brem_nurates[5:8] / eta_brem_WR[id_std][4:7]
ratio_nob = brem_nurates[5:8] / eta_brem_WR[id_nob][4:7]

for idx in range(3):
    # Eta coefficient
    axs[0][idx].plot(r, eta_brem_WR[id_std][idx+4], color = "tab:blue", ls = "-", label="WR - std")
    axs[0][idx].plot(r, eta_brem_WR[id_nob][idx+4], color = "tab:blue", ls = "--", label="WR - no bl")

    axs[0][idx].plot(r, brem_nurates[idx + 5], color = "tab:orange", label="NR")

    # Ratio
    axs[1][idx].plot(r, ratio_std[idx], color = "tab:blue", ls = "-" , label="WR - std")
    axs[1][idx].plot(r, ratio_nob[idx], color = "tab:blue", ls = "--", label="WR - no bl")

for ax in axs[0]:
    ax.set_yscale("log")
    ax.legend()

for ax in axs[:][1]:
    ax.axhline(1., color = "k", ls = "--", zorder=1)
    ax.set_xlabel(r"$r~[{\rm cm}]$")
    ax.set_yscale("log")

##########################################################

##################################
## Plot absoprtion coefficients ##
##################################

## Absoprtion on nucleons

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

axs[0][0].set_xscale("log")
plt.subplots_adjust(hspace=0.)

# Suptitle
plt.suptitle("Absorption coefficient for beta reactions")

# Title
axs[0][0].set_title(r"$\nu_{\rm e}$")
axs[0][1].set_title(r"$\bar{\nu}_{\rm e}$")
axs[0][2].set_title(r"$\nu_{\rm x}$")

# Ylabels
axs[0][0].set_ylabel(r"$\kappa_{\rm a}~[{\rm cm}^{-3}]$")
axs[1][0].set_ylabel("Ratio")

ratio_std = beta_nurates[11:13] / k_abs_WR[id_std][4:6]
ratio_nob = beta_nurates[11:13] / k_abs_WR[id_nob][4:6]

for idx in range(2):
    # Kappa coefficient
    axs[0][idx].plot(r, k_abs_WR[id_std][idx+4], color = "tab:blue", ls = "-" , label="WR - std")
    axs[0][idx].plot(r, k_abs_WR[id_nob][idx+4], color = "tab:blue", ls = "--", label="WR - no bl")

    axs[0][idx].plot(r, beta_nurates[idx+11], color = "tab:orange", label="NR")

    # Ratio
    axs[1][idx].plot(r, ratio_std[idx], color = "tab:blue", ls = "-" , label="WR - std")
    axs[1][idx].plot(r, ratio_nob[idx], color = "tab:blue", ls = "--", label="WR - no bl")

for ax in axs[0]:
    ax.set_yscale("log")
    ax.legend()

for ax in axs[:][1]:
    ax.axhline(1., color = "k", ls = "--", zorder=1)
    ax.set_xlabel(r"$r~[{\rm cm}]$")
    ax.set_yscale("log")

##########################################################

## Isoenergetic scattering

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

axs[0][0].set_xscale("log")
plt.subplots_adjust(hspace=0.)

# Suptitle
plt.suptitle("Absorption coefficient for isoenergetic scattering")

# Title
axs[0][0].set_title(r"$\nu_{\rm e}$")
axs[0][1].set_title(r"$\bar{\nu}_{\rm e}$")
axs[0][2].set_title(r"$\nu_{\rm x}$")

# Ylabels
axs[0][0].set_ylabel(r"$\kappa_{\rm s}~[{\rm cm}^{-3}]$")
axs[1][0].set_ylabel("Ratio")

ratio_std = iso_nurates[13:16] / k_sct_WR[id_std][4:7]
ratio_nob = iso_nurates[13:16] / k_sct_WR[id_nob][4:7]

for idx in range(3):
    # Kappa coefficient
    axs[0][idx].plot(r, k_sct_WR[id_std][idx+4], color = "tab:blue", ls = "-" , label="WR - std")
    axs[0][idx].plot(r, k_sct_WR[id_nob][idx+4], color = "tab:blue", ls = "--", label="WR - no bl")

    axs[0][idx].plot(r, iso_nurates[idx+13], color = "tab:orange", label="NR")

    # Ratio
    axs[1][idx].plot(r, ratio_std[idx], color = "tab:blue", ls = "-" , label="WR - std")
    axs[1][idx].plot(r, ratio_nob[idx], color = "tab:blue", ls = "--", label="WR - no bl")

for ax in axs[0]:
    ax.set_yscale("log")
    ax.legend()

for ax in axs[:][1]:
    ax.axhline(1., color = "k", ls = "--", zorder=1)
    ax.set_xlabel(r"$r~[{\rm cm}]$")
    ax.set_yscale("log")

plt.show()
