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

## Purpose of the library

BNS_NURATES can be exploited to compute both energy-dependent (spectral) and energy-integrated
(gray) emissivities and opacities for the most relevant neutrino-matter interactions in BNS
mergers. The definitions of such quantites can be found in [^fn1], which also discusses in
more details the design of the library.

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
and build
```
cmake --build . --target MWE
```
This will compile the minimal working example provided in the repo ```(mwe.cpp)```, which
evaluates spectral and gray neutrino rates for a single thermodynamic point. Thermodynamic
conditions extracted from a BNS merger simulation (used for the postprocessing
analysis in [^fn1]) are provided in the ```inputs/BNS``` folder for testing purposes.

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



## Python bindings

Bindings can be generated for calling BNS_NURATES functions directly from Python.
This is done through SWIG (tested version: 4.0.2). Once SWIG is installed on the system,
enter the ```bindings/``` folder
```
cd bindings
```
and run the following command
```
python3 setup.py build_ext --inplace
```
This creates a Python module named ```bnsnurates``` can that imported as
```
import bnsnurates as bns
```
To make the module available to your Python environment either add the ```bindings/``` folder
to your PYTHONPATH or install the module. ```bindings/test_bindings.py``` is a minimal working
example that can be used for testing the bindings.

## Included neutrino-matter interactions

BNS_NURATES implements the following neutrino reactions (see also [^fn1]):

- **Beta processes**

    Implemented as in [^fn2], [^fn3] under the assumption of zero-momentum transfer and non-relativistic nucleons. We account for the impact from
    relativistic mean-field effects as in [^fn3], [^fn4], and the correction due to the sum of phase-space, recoil and weak magnetism effects as in 
    [^fn5]. They include:

  * Neutrino absorptions on nucleons / $e^\pm$ captures
    ### $\nu_e + n \leftrightarrow e^- + p$
    ### $\bar{\nu}_e + p \leftrightarrow e^+ + n$

  * (Inverse) nucleon decay
    ### $\nu_e + n \leftrightarrow e^- + p$
    ### $\bar{\nu}_e + p \leftrightarrow e^+ + n$

- **Pair processes**

  * $e^+ e^-$ annihilation (implemented as in [^fn6])
    ### $e^+ + e^- \leftrightarrow \nu + \bar{\nu}$



  * Nucleon-nucleon bremsstrahlung
    - implementation 1: as [^fn7], including the in-medium modification from [^fn8]
    - implementation 2: as [^fn10]
    - implementation 3: as [^fn11]
    ### $N + N \leftrightarrow N + N + \nu + \bar{\nu}$

- **Scattering processes**

  * Isoenergetic scattering off nucleons (implemented as in [^fn2], including phase-space, recoil and weak magnetism effects as in 
    [^fn5])
    ### $\nu + N \rightarrow \nu + N$

  * Inelastic scattering off $e^\pm$ (implemented as in [^fn2], [^fn9])
    ### $\nu + e^\pm \rightarrow \nu + e^\pm$


[^fn1]: [L. Chiesa et al., Phys. Rev. D 111, 063053 (2025)](https://doi.org/10.1103/PhysRevD.111.063053)
[^fn2]: [S. W. Bruenn, Astrophys. J. Suppl. Ser. 58, 771 (1985)](https://doi.org/10.1086/191056)
[^fn3]: [M. Oertel, A. Pascal, M. Mancini, and J. Novak, Phys. Rev. C 102, 035802 (2020)](https://doi.org/10.1103/PhysRevC.102.035802)
[^fn4]: [M. Hempel, Phys. Rev. C 91, 055807 (2015)](https://doi.org/10.1103/PhysRevC.91.055807)
[^fn5]: [C. J. Horowitz, Phys. Rev. D 65, 043001 (2002)](https://doi.org/10.1103/PhysRevD.65.043001)
[^fn6]: [J. A. Pons, J. A. Miralles, and J. M. Ibanez, Astron. Astrophys. Suppl. Ser. 129, 343 (1998)](https://doi.org/10.1051/aas:1998189)
[^fn7]: [S. Hannestad and G. Raffelt, Astrophys. J. 507, 339 (1998)](https://doi.org/10.1086/306303)
[^fn8]: [T. Fischer, Astron. Astrophys. 593, A103 (2016)](https://doi.org/10.1051/0004-6361/201628991)
[^fn9]: [A. Mezzacappa and S. W. Bruenn, Astrophys. J. 410, 740 (1993)](https://doi.org/10.1086/172791)
[^fn10]:[A. Burrows, S. Reddy, and T. A. Thompson, Nucl. Phys. A777, 356 (2006)](https://doi.org/10.1016/j.nuclphysa.2004.06.012)
[^fn11]:[Guo G. and Martínez-Pinedo G., ApJ 887 58 (2019)](https://doi.org/10.3847/1538-4357/ab536d)
