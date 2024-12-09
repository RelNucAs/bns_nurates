# BNS_NURATES: A GPU friendly library for neutrino opacities written in c++ with GPU support

BNS_NURATES requires the following dependencies:
- A C++ compiler which supports c++17
- The GNU scientific library
- Kokkos
- OpenMP

## Installation

BNS_NURATES is a header-only library, so it can be simply included into your
code, provided the requirements listed above are satisfied. The various headers
are in the `include` directory.


## Stand alone build and testing

Clone the repository recursively to also clone Kokkos which is available as a submodule:
```
git clone --recursive https://github.com/RelNucAs/bns_nurates.git
cd bns_nurates
```
Create a build directory:
```
mkdir build
cd build
```
By default Kokkos and OpenMP are enabled. To compile for a CPU based system run
```
cmake ../
```
and build nurates
```
cmake --build . --target nurates
```
If the code is intended for Nvidia GPUs (say, A100 in this case), enable CUDA during compilation
```
cmake -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_AMPERE80=ON ../
```

Similarly, to disable OpenMP, add the flag ```-DENABLE_OPENMP=OFF``` to cmake.

## Using bns_nurates as a thorn for the Einsten Toolkit
bns_nurates is also meant for use as a thorn ```Weakrates2``` for Cactus with
THC. To do this, one must export these files:

```
python utils/export_thc.py /path/to/bns_nurates/ /path/to/destination/
```

Here ```/path/to/bns_nurates/``` is the full path for the top-level bns_nurates
directory and ```/path/to/destination/``` is the full path to the ```src```
folder insider ```Weakrates2```


##### A C library for neutrino opacities in the context of binary neutron star mergers.

These opacities are based on Burrows (2004)[^fn1], Bruenn (1985)[^fn2] and
references therein.

The following reactions are implemented:
- Neutrino absorption on protons Eq:C13 from Bruenn (1985) with blocking factor
  and weak magnetism from Burrows (2004) included.
- Neutrino absorption on neutrons Eq:C19 from Bruenn (1985) with blocking factor
  and weak magnetism from Burrows (2004) included.
- Pair production and annihilation of neutrinos Eq:2, Eq:3 and emissivity, mean
  free path from Eq:4 of Pons et. al. (1998)[^fn3] are included.

[^fn1]: [Burrows, Reddy and Thompson, Neutrino opacities in nuclear matter,
    Nucl.Phys. A 777 356-394
    (2006)](https://doi.org/10.1016/j.nuclphysa.2004.06.012)
[^fn2]: [Bruenn, Stellar core collapse - Numerical model and infall epoch,
    Astrophysical Journal Supplement Series 58 771-841
    (1985)](https://doi.org/10.1086/191056)
[^fn3]: [Pons, Miralles and Ibáñez, Legendre expansion of the kernel nu nubar ->
e+e-: Influence of high order terms, Astron. Astrophys. Suppl. Ser. 129, 343-351
(1998)](https://doi.org/10.1051/aas:1998189) TODO: update this list
