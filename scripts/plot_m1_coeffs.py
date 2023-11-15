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
    

## Read WeakRates data
def read_WeakRates_data(filename):
    foldername = "../inputs/wr_standard/"
    data_1 = np.loadtxt(foldername + filename, comments="#", delimiter=" ", unpack=True)

    foldername = "../inputs/wr_noblock/"
    data_2 = np.loadtxt(foldername + filename, comments="#", delimiter=" ", unpack=True)

    return [data_1, data_2]

eta_pair_WR = read_WeakRates_data("eta_pair.dat")
eta_brem_WR = read_WeakRates_data("eta_Bremsstrahlung.dat")
eta_plas_WR = read_WeakRates_data("eta_plasmon.dat")
eta_beta_WR = read_WeakRates_data("eta_betadecay.dat")

id_std = 0 # index for standard WeakRates data
id_nob = 1 # index for modified (no blocking) WeakRates data

r = eta_beta_WR[id_std][0] # radius of the profile


## Read bns_nurates data
def read_nurates_data(filename):
    foldername = "../tests/tests_opacities_m1/output/"
    if args.mod: filename = filename.replace(".txt","_mod.txt")
    data = np.loadtxt(foldername + filename, comments="#", delimiter=" ", unpack=True)
    return data

k_abs_WR = read_WeakRates_data("k_absorption_np.dat")
k_sct_WR = read_WeakRates_data("k_scattering_np.dat")

beta_nurates = read_nurates_data("m1_opacities_abs_em.txt")
iso_nurates  = read_nurates_data("m1_opacities_isoscatt.txt")
pair_nurates = read_nurates_data("m1_opacities_pair.txt")
brem_nurates = read_nurates_data("m1_opacities_brem.txt")

## Generic function for data plotting
def make_plot(axs, WR_data, nurates_data, id_type, id_coeff):
    # Define index for reading data correctly
    if not (0 <= id_type <= 1):
        print("id_type must be either 0 (number) or 1 (energy)")
        exit()

    if not (0 <= id_coeff <= 2):
        print("id_coef must be either 0 (eta), 1 (kappa_a) or 2 (kappa_s)")
        exit()

    id_WR = 3 * id_type + 1 
    id_NR = 6 * id_coeff + 3 * id_type + 2 # two additional columns for r and diff_distr

    if (id_coeff == 2):
        id_NR = 6 * id_coeff + 2 - 1 # minus one since diff_distr column is missing for iso scattering

    axs[0][0].set_xscale("log")
    plt.subplots_adjust(hspace=0.)
    
    # Title
    axs[0][0].set_title(r"$\nu_{\rm e}$")
    axs[0][1].set_title(r"$\bar{\nu}_{\rm e}$")
    axs[0][2].set_title(r"$\nu_{\rm x}$")

    # Ylabel
    axs[1][0].set_ylabel("Ratio")

    # Compute ratio
    ratio_std = nurates_data[id_NR:id_NR+3] / WR_data[id_std][id_WR:id_WR+3]
    ratio_nob = nurates_data[id_NR:id_NR+3] / WR_data[id_nob][id_WR:id_WR+3]
    
    # Plot data
    for idx in range(3):
        # Rates
        axs[0][idx].plot(r, WR_data[id_std][id_WR+idx], color = "tab:blue", ls = "-", label="WR - std")
        axs[0][idx].plot(r, WR_data[id_nob][id_WR+idx], color = "tab:blue", ls = "--", label="WR - no bl")

        axs[0][idx].plot(r, nurates_data[id_NR+idx], color = "tab:orange", label="NR")

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

################################
## Plot emission coefficients ##
################################

## id_type = 0 (number), 1 (energy)
## id_coeff = 0 (eta), 1 (kappa_a), 2 (kappa_s)

## Beta decay - number

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

# Suptitle
plt.suptitle("Emissivity coefficient for beta reactions - number")

# Ylabels
axs[0][0].set_ylabel(r"$\eta^0~[{\rm MeV}\,{\rm cm}^{-3}\,{\rm s}^{-1}]$")

# Make plot
make_plot(axs, eta_beta_WR, beta_nurates, 0, 0)

##########################################################

## Beta decay - energy

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

# Suptitle
plt.suptitle("Emissivity coefficient for beta reactions - energy")

# Ylabels
axs[0][0].set_ylabel(r"$\eta~[{\rm MeV}\,{\rm cm}^{-3}\,{\rm s}^{-1}]$")

# Make plot
make_plot(axs, eta_beta_WR, beta_nurates, 1, 0)


##########################################################

## Pair process - number 

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

# Suptitle
plt.suptitle(r"Emissivity coefficient for $e^\pm$ annihilation - number")

# Ylabels
axs[0][0].set_ylabel(r"$\eta^0~[{\rm MeV}\,{\rm cm}^{-3}\,{\rm s}^{-1}]$")

# Make plot
make_plot(axs, eta_pair_WR, pair_nurates, 0, 0)

##########################################################

## Pair process - energy 

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

# Suptitle
plt.suptitle(r"Emissivity coefficient for $e^\pm$ annihilation - energy")

# Ylabels
axs[0][0].set_ylabel(r"$\eta~[{\rm MeV}\,{\rm cm}^{-3}\,{\rm s}^{-1}]$")

# Make plot
make_plot(axs, eta_pair_WR, pair_nurates, 1, 0)

##########################################################

## Bremsstrahlung process - number

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

# Suptitle
plt.suptitle(r"Emissivity coefficient for NN bremsstrahlung - number")

# Ylabels
axs[0][0].set_ylabel(r"$\eta^0~[{\rm MeV}\,{\rm cm}^{-3}\,{\rm s}^{-1}]$")

# Make plot
make_plot(axs, eta_brem_WR, brem_nurates, 0, 0)

##########################################################

## Bremsstrahlung process - energy

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

# Suptitle
plt.suptitle(r"Emissivity coefficient for NN bremsstrahlung - energy")

# Ylabels
axs[0][0].set_ylabel(r"$\eta~[{\rm MeV}\,{\rm cm}^{-3}\,{\rm s}^{-1}]$")

# Make plot
make_plot(axs, eta_brem_WR, brem_nurates, 1, 0)

##########################################################

##################################
## Plot absoprtion coefficients ##
##################################

## Absoprtion on nucleons - number

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

# Suptitle
plt.suptitle(r"Absorption coefficient for beta reactions - number")

# Ylabels
axs[0][0].set_ylabel(r"$\kappa^0_{\rm a}~[{\rm MeV}\,{\rm cm}^{-3}\,{\rm s}^{-1}]$")

# Make plot
make_plot(axs, k_abs_WR, beta_nurates, 0, 1)

##########################################################

## Absoprtion on nucleons - energy

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

# Suptitle
plt.suptitle(r"Absorption coefficient for beta reactions - energy")

# Ylabels
axs[0][0].set_ylabel(r"$\kappa_{\rm a}~[{\rm MeV}\,{\rm cm}^{-3}\,{\rm s}^{-1}]$")

# Make plot
make_plot(axs, k_abs_WR, beta_nurates, 1, 1)

##########################################################

## Isoenergetic scattering - energy

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]})

# Suptitle
plt.suptitle("Absorption coefficient for isoenergetic scattering - energy")

# Ylabels
axs[0][0].set_ylabel(r"$\kappa_{\rm s}~[{\rm MeV}\,{\rm cm}^{-3}\,{\rm s}^{-1}]$")

# Make plot
make_plot(axs, k_sct_WR, iso_nurates, 1, 2)

plt.show()