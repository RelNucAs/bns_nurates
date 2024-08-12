# bns_nurates: A performance portable library for neutrino opacities in Kokkos and C++

By default, bns_nurates compiles with OpenMP and Kokkos. To prevent this, add flags -DENABLE_KOKKOS=OFF and -DENABLE_OPENMP=OFF to cmake.

##### A C library for neutrino opacities in the context of binary neutron star mergers.

These opacities are based on Burrows (2004)[^fn1], Bruenn (1985)[^fn2] and references therein.

The following reactions are implemented:
- Neutrino absorption on protons Eq:C13 from Bruenn (1985) with blocking factor and weak magnetism from Burrows (2004) included.
- Neutrino absorption on neutrons Eq:C19 from Bruenn (1985) with blocking factor and weak magnetism from Burrows (2004) included.
- Pair production and annihilation of neutrinos Eq:2, Eq:3 and emissivity, mean free path from Eq:4 of Pons et. al. (1998)[^fn3] are included.

[^fn1]: [Burrows, Reddy and Thompson, Neutrino opacities in nuclear matter, Nucl.Phys. A 777 356-394 (2006)](https://doi.org/10.1016/j.nuclphysa.2004.06.012)
[^fn2]: [Bruenn, Stellar core collapse - Numerical model and infall epoch, Astrophysical Journal Supplement Series 58 771-841 (1985)](https://doi.org/10.1086/191056)
[^fn3]: [Pons, Miralles and Ibáñez, Legendre expansion of the kernel nu nubar -> e+e-: Influence of high order terms, Astron. Astrophys. Suppl. Ser. 129, 343-351 (1998)](https://doi.org/10.1051/aas:1998189)
