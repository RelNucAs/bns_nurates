# NuO-BNSlib

##### A C library for neutrino opacities in the context of binary neutron star mergers.

These opacities are based on Burrows (2004)[^fn1], Bruenn (1985)[^fn2] and references therein.

The following reactions are implemented:
- Neutrino absorption on protons Eq:C13 from Bruenn (1985) with blocking factor and weak magnetism from Burrows (2004) included.
- Neutrino absorption on neutrons Eq:C19 from Bruenn (1985) with blocking factor and weak magnetism from Burrows (2004) included.
- Pair production and annihilation of neutrinos Eq:2, Eq:3 and emissivity, mean free path from Eq:4 of Pons et. al. (1998)[^fn3] are included.

Code description:

For the moment main.cpp computes the (energy-integrated) emissivity and number emissivity for (anti)electron neutrinos along a 1D CCSN profile.

Use Makefile to compile the main file.

[^fn1]: Burrows, Reddy and Thompson, Neutrino opacities in nuclear matter, Nucl.Phys. A777 356-394 (2006)
[^fn2]: Bruenn, Stellar core collapse - Numerical model and infall epoch
[^fn3]: Pons, Miralles and Ibáñez, Legendre expansion of the kernel: Influence of high order terms, Astron. Astrophys. Suppl. Ser. 129, 343-351 (1998)
