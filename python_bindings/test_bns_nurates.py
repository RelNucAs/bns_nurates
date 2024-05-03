import numpy as np
import bns_nurates as bns

#######################################
# Testing bns_nurates python bindings #
#######################################

# Define physical constants
c = 2.997924562E10
MeV = 1.602176634E-6

print("Testing python bindings generated with SWIG for the bns_nurates library...")

# Define and print thermodynamics conditions
m_ref =  931.494061
rho  = 1435091918127104.0   # g cm-3
temp = 51.150005340576172   # MeV
ye   = 0.21549025177955627  
xn   = 0.78450974821565389
xp   = 0.21549025178434600
mu_e = 324.00790839794519   # MeV
mu_n = 585.93269122709501   # MeV
mu_p = 419.35214965870006   # MeV

nb = rho * c * c / (m_ref * MeV)

print("Thermodynamics conditions:")
print("rho  = %lf g cm-3" %rho)
print("nb   = %lf cm-3  " %nb)
print("T    = %lf MeV   " %temp)
print("Ye   = %lf       " %ye)
print("Yn   = %lf       " %xn)
print("Yp   = %lf       " %xp)
print("mu_e = %lf MeV   " %mu_e)
print("mu_n = %lf MeV   " %mu_n)
print("mu_p = %lf MeV   " %mu_p)
print("")

# Populate EOS structure
eos_pars = bns.MyEOSParams()
eos_pars.nb = nb
eos_pars.temp = temp
eos_pars.ye = ye
eos_pars.yn = xn
eos_pars.yp = xp
eos_pars.mu_e = mu_e
eos_pars.mu_n = mu_n
eos_pars.mu_p = mu_p

# Set reactions
opacity_flags = {'use_abs_em': True, 'use_pair': True, 'use_brem': False, 'use_inelastic_scatt': True, 'use_iso': True}
bns.print_reactions(opacity_flags)

# Set reactions
opacity_pars = {'use_dU': False, 'use_dm_eff': False, 'use_WM_ab': False, 'use_WM_sc': False, 'use_decay': False, 'neglect_blocking': False}
bns.print_corrections(opacity_pars)

# Set distribution function and radiation quantities
print("Neutrino distribution function reconstruction: 'equilibrium'\n")
distr_pars = bns.NuEquilibriumParams(eos_pars)

m1_pars = bns.M1Quantities()
Jvals = bns.NuEnergy(distr_pars)
nvals = bns.NuNumber(distr_pars)
m1_pars.n = nvals.integrand[:4]
m1_pars.J = Jvals.integrand[:4]

# Populate global structure
grey_pars = bns.GreyOpacityParams()
grey_pars.eos_pars = eos_pars
grey_pars.opacity_flags = opacity_flags
grey_pars.opacity_pars = opacity_pars
grey_pars.distr_pars = distr_pars
grey_pars.m1_pars = m1_pars

# Compute quadrature
quad = bns.cvar.quadrature_default
bns.GaussLegendreMultiD(quad)
print("Number of quadrature points: %d\n" %quad.nx)

# Set neutrino energy
e_nu = 0.5 # MeV
print("Neutrino energy: %lf MeV\n" %e_nu)

# Compute non stimulated spectral opacities
spec_opacities = bns.ComputeSpectralOpacitiesNotStimulatedAbs(e_nu, quad, grey_pars)
print("Spectral rates -> non stimulated")
bns.print_spectral_rates(spec_opacities)

# Compute stimulated spectral opacities
spec_opacities_stim = bns.ComputeSpectralOpacitiesStimulatedAbs(e_nu, quad, grey_pars)
print("Spectral rates -> stimulated")
bns.print_spectral_rates(spec_opacities_stim)

# Compute non stimulated integrated opacities
int_opacities = bns.ComputeM1OpacitiesNotStimulated(quad, quad, grey_pars)
print("Integrated rates -> non stimulated")
bns.print_integrated_rates(int_opacities)

# Compute stimulated integrated opacities
int_opacities_stim = bns.ComputeM1Opacities(quad, quad, grey_pars)
print("Integrated rates -> stimulated")
bns.print_integrated_rates(int_opacities_stim)
