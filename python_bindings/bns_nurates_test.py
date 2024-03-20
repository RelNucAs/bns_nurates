import bns_nurates

c = 2.997924562E10
MeV = 1.602176634E-6

rho = 1435091918127104.0
temp =  51.150005340576172
ye = 0.21549025177955627
xn = 0.78450974821565389
xp = 0.21549025178434600
abar = 1.0000000000000000
zbar = 1.0000000000000000
mu_e = 324.00790839794519
mu_n = 585.93269122709501
mu_p = 419.35214965870006
mu_hat = 166.58054156839492

m_ref =  931.494061

n_leg = 50

opacity_flags = bns_nurates.OpacityFlags()
opacity_flags.use_abs_em = 0
opacity_flags.use_iso = 0
opacity_flags.use_brem = 1
opacity_flags.use_pair = 0
opacity_flags.use_inelastic_scatt = 0

nb = rho * c * c / (m_ref * MeV)

eos_pars = bns_nurates.MyEOSParams()

eos_pars.nb = nb
eos_pars.temp = temp
eos_pars.ye = ye
eos_pars.yn = xn
eos_pars.yp = xp
eos_pars.mu_e = mu_e
eos_pars.mu_n = mu_n
eos_pars.mu_p = mu_p

e_nu = 5.00000000e-01

quad_1d = bns_nurates.MyQuadrature()
quad_1d.nx = n_leg
quad_1d.dim = 1
quad_1d.type = bns_nurates.kGauleg
quad_1d.x1 = 0
quad_1d.x2 = 1

bns_nurates.GaussLegendreMultiD(quad_1d)

distr_pars = bns_nurates.NuEquilibriumParams(eos_pars)

m1_pars = bns_nurates.M1Quantities()

Jvals = bns_nurates.NuEnergy(distr_pars)
nvals = bns_nurates.NuNumber(distr_pars)

m1_pars.n = nvals.integrand
m1_pars.J = Jvals.integrand

grey_pars = bns_nurates.GreyOpacityParams()
grey_pars.eos_pars = eos_pars
grey_pars.opacity_flags = opacity_flags
grey_pars.opacity_params = bns_nurates.cvar.opacity_params_default_none
grey_pars.distr_pars = distr_pars
grey_pars.m1_pars = m1_pars


spec_opacities = bns_nurates.ComputeSpectralOpacitiesNotStimulatedAbs(e_nu, quad_1d, grey_pars)
m1_opacities = bns_nurates.ComputeM1Opacities(quad_1d, quad_1d, grey_pars)
print(m1_opacities.eta)

print(dir(m1_opacities))
print(spec_opacities.j)
print(spec_opacities.j_s)
print(spec_opacities.kappa)
print(spec_opacities.kappa_s)
spec_opacities.j = [0.,1.,2.,3.]
print(spec_opacities.j)
