import bnsnurates as bns

#########################################
# Test of BNS_NURATES's Python bindings #
#########################################

## Fix the neutrino energy for computation of spectral rates
nu_energy = 10. # [MeV]

## Input thermodynamic quantities (corresponding to point A in Chiesa+25 PRD)
##  N.B.: chemical potentials include the rest mass contribution
nb   = 4.208366627847035e+38  # Baryon number density [cm-3]           
temp = 12.406403541564941     # Temperature [MeV]
ye   = 0.07158458232879639    # Electron fraction
xn   = 1. - ye                # Neutron fraction
xp   = ye                     # Proton fraction
mu_e = 187.1814489            # Electron chemical potential [MeV]
mu_p = 1011.01797737          # Proton chemical potential [MeV]
mu_n = 1221.59013681          # Neutron chemical potential [MeV]
dU = 18.92714728              # Nucleon interaction potential difference (Un-Up) [MeV]
mp_eff = 278.87162217         # Proton effective mass [MeV]
mn_eff = 280.16495513         # Neutron effective mass [MeV]
dm = mn_eff - mp_eff          # Nucleon effective mass difference [MeV]

## Input gray neutrino quantities (library supports 4 neutrino species)
## Defined as [nue, anue, nux, anux]
## N.B.: These will be used in PART 2 for reconstructing the neutrino distribution functions
##       and as normalization factors for energy-averaged opacities
n_m1 = [3.739749408027436e+33, 1.2174961961689319e+35, 2.2438496448164613e+34, 2.2438496448164613e+34]  # Neutrino number densities [cm-3]
J_m1 = [1.246583136009145e+35,  5.360307484839323e+36,  8.726081952064015e+35,  8.726081952064015e+35]  # Neutrino energy densities [MeV cm-3]
chi_m1  = [1./3., 1./3., 1./3., 1./3.]  # Eddington factor

## Pass thermodynamic conditions
eos_pars = bns.MyEOSParams()
eos_pars.nb     = nb * 1e-21  # Set baryon number density (in baryon/nm^3)
eos_pars.temp   = temp        # Set temperature (in MeV)
eos_pars.ye     = ye          # Set electron fraction
eos_pars.yn     = xn          # Set neutron fraction
eos_pars.yp     = xp          # Set proton fraction
eos_pars.mu_e   = mu_e        # Set electron chemical potential (in MeV)
eos_pars.mu_n   = mu_n        # Set neutron  chemical potential (in MeV)
eos_pars.mu_p   = mu_p        # Set proton   chemical potential (in MeV) (NOTE: reactions only depend on the difference between mu_n and mu_p)
eos_pars.dU     = dU          # Set effective potential difference (in MeV)
eos_pars.dm_eff = dm          # Set effective mass difference (in MeV)

## Create a quadrature struct and populate it with data relative to
## Gauss-Legendre quadrature
quad = bns.cvar.quadrature_default
quad.nx = 6  # Number of quadrature points
bns.GaussLegendre(quad)
print("Number of quadrature points: %d\n" %quad.nx)

## Select and print active reactions
'''
- use_abs_em          : neutrino abs on nucleons - e+/e- captures
- use_brem            : nucleon-nucleon bremsstrahlung
- use_pair            : e+ e- annihilation (and inverse)
- use_iso             : isoenergetic scattering off nucleons
- use_inelastic_scatt : inelastic scattering off e+/e-
'''
opacity_flags = {'use_abs_em': True, 'use_pair': True, 'use_brem': True, 'use_inelastic_scatt': True, 'use_iso': True}
bns.print_reactions(opacity_flags)

## Select and print active corrections
'''
- use_dU             : dU correction to beta processes
- use_dm_eff         : effective dm correction to beta processes
- use_WM_ab          : weak magnetism correction to beta processes
- use_WM_sc          : weak magnetism correction to scattering off nucleons
- use_decay          : include (inverse) nucleon decays in beta processes
- use_BRT_brem       : implement bremsstrahlung as in Burrows+2006
- neglect_blocking   : neglect blocking factors
- use_NN_medium_corr : medium correction to bremsstrahlung as Fischer+16
'''
opacity_pars = {'use_dU': True, 'use_dm_eff': False, 'use_WM_ab': True, 'use_WM_sc': True, 'use_decay': True, 'use_BRT_brem': False, 'neglect_blocking': False, 'use_NN_medium_corr': True}
bns.print_corrections(opacity_pars)

print("\nInput parameters")
print("----------------")
print("Neutrino energy (for spectral rates)               : %13.6e (MeV)" %nu_energy)
print("Baryon number density                              : %13.6e (cm^-3)" %nb)
print("Temperature                                        : %13.6e (MeV)" %temp)
print("Electron fraction                                  : %13.3f" %ye)
print("Relativistic electron chemical potential           : %13.6e (MeV)" %mu_e)
print("Relativistic proton chemical potential             : %13.6e (MeV)" %mu_p)
print("Relativistic neutron chemical potential            : %13.6e (MeV)" %mu_n)
print("Effective nucleon mass difference (n - p)          : %13.6e (MeV)" %dm)
print("Effective nucleon potential difference (n - p)     : %13.6e (MeV)" %dU)
print("Electron neutrinos number density 'n'              : %13.6e (cm^-3)" %n_m1[0])
print("Electron antineutrinos number density 'n'          : %13.6e (cm^-3)" %n_m1[1])
print("Heavy-type neutrinos number density 'n'            : %13.6e (cm^-3)" %n_m1[2])
print("Heavy-type antineutrinos number density 'n'        : %13.6e (cm^-3)" %n_m1[3])
print("Electron neutrinos energy density 'J'              : %13.6e (MeV cm^-3)" %J_m1[0])
print("Electron antineutrinos energy density 'J'          : %13.6e (MeV cm^-3)" %J_m1[1])
print("Heavy-type neutrinos energy density 'J'            : %13.6e (MeV cm^-3)" %J_m1[2])
print("Heavy-type antineutrinos energy density 'J'        : %13.6e (MeV cm^-3)" %J_m1[3])
print("Electron neutrinos Eddington parameter 'chi'       : %13.11f" %chi_m1[0])
print("Electron antineutrinos Eddington parameter 'chi'   : %13.11f" %chi_m1[1])
print("Heavy-type neutrinos Eddington parameter 'chi'     : %13.11f" %chi_m1[2])
print("Heavy-type antineutrinos Eddington parameter 'chi' : %13.11f" %chi_m1[3])
print("\n")

############################################################
## PART 1: compute rates assuming equilibrium with matter ## 
############################################################

## Compute neutrino distribution parameters assuming equilibrium
distr_pars = bns.NuEquilibriumParams(eos_pars)

## Compute gray neutrino number and energy densities assuming equilibrium
## N.B.: they are required for the normalization factor of energy-averaged opacities,
##       they do not enter in the computation of spectral rates
m1_pars = bns.M1Quantities()
bns.ComputeM1DensitiesEq(eos_pars, distr_pars, m1_pars)
m1_pars.chi = [1./3., 1./3., 1./3., 1./3.]
m1_pars.J = [x * bns.kBS_MeV for x in m1_pars.J] # restore to input units (from MeV nm^-3 to g s^-2 nm^-1)

## Restore units of cm^-3 and MeV cm^-3 for printing
n_eq = [x * 1e21 for x in m1_pars.n]
J_eq = [x * 1e21 for x in m1_pars.J]
chi_eq = m1_pars.chi

print("Reconstructed neutrino densities assuming equilibrium");
print("-----------------------------------------------------");
print("    nue            anue            nux            anux");
print("n  %13.6e  %13.6e   %13.6e  %13.6e     (cm^-3)" %(n_eq[0], n_eq[1], n_eq[2], n_eq[3]))
print("J  %13.6e  %13.6e   %13.6e  %13.6e (MeV cm^-3)" %(J_eq[0], J_eq[1], J_eq[2], J_eq[3]))
print("chi %12.10f   %12.10f    %12.10f   %12.10f"   %(chi_eq[0], chi_eq[1], chi_eq[2], chi_eq[3]))
print("")

## Populate global structure
grey_pars = bns.GreyOpacityParams()
grey_pars.eos_pars = eos_pars
grey_pars.opacity_flags = opacity_flags
grey_pars.opacity_pars = opacity_pars
grey_pars.distr_pars = distr_pars
grey_pars.m1_pars = m1_pars

## Compute and output spectral emissivities and inverse mean free paths (not in the stimulated
## absorption formalism)
spectral_rates = bns.ComputeSpectralOpacitiesNotStimulatedAbs(nu_energy, quad, grey_pars) #bns.ComputeSpectralOpacitiesStimulatedAbs() for stimulated absorption 

## The numerical factors restore usual units (see output)
spectral_rates['kappa']   = [x * 1e7 for x in spectral_rates['kappa']]
spectral_rates['kappa_s'] = [x * 1e7 for x in spectral_rates['kappa_s']]

## N.B.: 'j' and 'kappa' are the spectral emissivities and inverse mean free paths for the
##       sum of the active inelastic reactions, 'j_s' and 'kappa_s' are the equivalent for
##       the elastic scattering off nucleons
print("Spectral rates assuming equilibrium")
print("------------------------------")
bns.print_spectral_rates(spectral_rates)

## Compute and output gray emissivities and opacities (Eqs. (19)-(23) in Chiesa+25 PRD)
gray_rates = bns.ComputeM1Opacities(quad, quad, grey_pars)

## The numerical factors restore usual units (see output)
gray_rates['eta']       = [x * 1e21 for x in gray_rates['eta']]
gray_rates['eta_0']     = [x * 1e21 for x in gray_rates['eta_0']]
gray_rates['kappa_a']   = [x * 1e7  for x in gray_rates['kappa_a']]
gray_rates['kappa_0_a'] = [x * 1e7  for x in gray_rates['kappa_0_a']]
gray_rates['kappa_0_s'] = [x * 1e7  for x in gray_rates['kappa_s']]

print("Gray rates assuming equilibrium")
print("------------------------------")
bns.print_integrated_rates(gray_rates)

####################################################################
## PART 2: compute rates reconstructing the neutrino distribution ##
##         function from neutrino number and energy densities     ##
####################################################################

## Set neutrino number and energy densities
m1_pars.n = [x * 1e-21 for x in n_m1]               # nm^-3
m1_pars.J = [x * 1e-21 * bns.kBS_MeV for x in J_m1] # g s^-2 nm^-1
m1_pars.chi = chi_m1

## Compute neutrino distribution parameters from neutrino number/energy densities
distr_pars = bns.CalculateDistrParamsFromM1(m1_pars, eos_pars)

## Assign neutrino input to global structure
grey_pars.m1_pars    = m1_pars
grey_pars.distr_pars = distr_pars

## Compute and output spectral emissivities and inverse mean free paths (not in the stimulated
## absorption formalism)
spectral_rates = bns.ComputeSpectralOpacitiesNotStimulatedAbs(nu_energy, quad, grey_pars)

## The numerical factors restore usual units (see output)
spectral_rates['kappa']   = [x * 1e7 for x in spectral_rates['kappa']]
spectral_rates['kappa_s'] = [x * 1e7 for x in spectral_rates['kappa_s']]

print("Spectral rates reconstructing distribution function");
print("------------------------------");
bns.print_spectral_rates(spectral_rates)

## Compute and output gray emissivities and opacities (Eqs. (19)-(23) in Chiesa+25 PRD)
gray_rates = bns.ComputeM1Opacities(quad, quad, grey_pars)

## The numerical factors restore usual units (see output)
gray_rates['eta']       = [x * 1e21 for x in gray_rates['eta']]
gray_rates['eta_0']     = [x * 1e21 for x in gray_rates['eta_0']]
gray_rates['kappa_a']   = [x * 1e7  for x in gray_rates['kappa_a']]
gray_rates['kappa_0_a'] = [x * 1e7  for x in gray_rates['kappa_0_a']]
gray_rates['kappa_0_s'] = [x * 1e7  for x in gray_rates['kappa_s']]

print("Gray rates reconstructing distribution function");
print("----------------------------------------------");
bns.print_integrated_rates(gray_rates)

print("Units\n"
      "-----\n"
      "Spectral emissivity 'j'/'j_s'   :           s^-1\n"
      "Spectral imfp 'kappa'/'kappa_s' :          cm^-1\n"
      "Gray number emissivity 'eta0'   :     cm^-3 s^-1\n"
      "Gray energy emissivity 'eta1'   : MeV cm^-3 s^-1\n"
      "Gray number opacity 'kappa0'    :     cm^-1 s^-1\n"
      "Gray energy opacity 'kappa1'    : MeV cm^-1 s^-1\n"
      "Gray scattering opacity 'scat1' : MeV cm^-1 s^-1\n")


