// bnsnurates.i

// @TODO: integrate with numpy.i

%module bnsnurates

%{
#include "bns_nurates.hpp"
#include "constants.hpp"
#include "distribution.hpp"
#include "eos.hpp"
#include "formulas_beta_iso.hpp"
#include "functions.hpp"
#include "integration.hpp"
#include "kernel_brem.hpp"
#include "kernel_nes.hpp"
#include "kernel_pair.hpp"
#include "kernels.hpp"
#include "m1_opacities.hpp"
#include "weak_magnetism.hpp"
%}

// Define __attribute__ directive away to nothing
#define __attribute__(x) 

// Exclude the following typedef struct from wrapping
%ignore M1MatrixKokkos2D;

///////////////////////////
// Define input typemaps //
///////////////////////////

// Python list to BS_REAL array
%typemap(in) BS_REAL [ANY] (BS_REAL tmp[$1_dim0]) {
    int i;
    if (!PySequence_Check($input)) {
            PyErr_SetString(PyExc_ValueError, "Expected a sequence");
            SWIG_fail;
    }
    if (PySequence_Length($input) != $1_dim0) {
            PyErr_SetString(PyExc_ValueError, "Size mismatch. Expected $1_dim0 elements");
            SWIG_fail;
    }
    for (i = 0; i < $1_dim0; i++) {
            PyObject *o = PySequence_GetItem($input, i);
            if (PyNumber_Check(o)) {
                    tmp[i] = (double) PyFloat_AsDouble(o);
            } else {
                    PyErr_SetString(PyExc_ValueError, "Sequence elements must be numbers");
                    SWIG_fail;
            }
    }
    $1 = tmp;
}

// Python dict to C++ OpacityFlags struct
%typemap(in) OpacityFlags * (OpacityFlags tmp) {
    tmp.use_abs_em = (int) PyInt_AsLong(PyDict_GetItemString($input, "use_abs_em"));
    tmp.use_pair = (int) PyInt_AsLong(PyDict_GetItemString($input, "use_pair"));
    tmp.use_brem = (int) PyInt_AsLong(PyDict_GetItemString($input, "use_brem"));
    tmp.use_inelastic_scatt = (int) PyInt_AsLong(PyDict_GetItemString($input, "use_inelastic_scatt"));
    tmp.use_iso = (int) PyInt_AsLong(PyDict_GetItemString($input, "use_iso"));
    
    $1 = &tmp;
}

// Python dict to C++ OpacityParams struct
%typemap(in) OpacityParams * (OpacityParams tmp) {
    tmp.use_dU = PyInt_AsLong(PyDict_GetItemString($input, "use_dU"));
    tmp.use_dm_eff = PyInt_AsLong(PyDict_GetItemString($input, "use_dm_eff"));
    tmp.use_WM_ab = PyInt_AsLong(PyDict_GetItemString($input, "use_WM_ab"));
    tmp.use_WM_sc = PyInt_AsLong(PyDict_GetItemString($input, "use_WM_sc"));
    tmp.use_decay = PyInt_AsLong(PyDict_GetItemString($input, "use_decay"));
    tmp.use_BRT_brem = PyInt_AsLong(PyDict_GetItemString($input, "use_BRT_brem"));
    tmp.use_NN_medium_corr = PyInt_AsLong(PyDict_GetItemString($input, "use_NN_medium_corr"));
    tmp.neglect_blocking = PyInt_AsLong(PyDict_GetItemString($input, "neglect_blocking"));

    $1 = &tmp;
}

////////////////////////////
// Define output typemaps //
////////////////////////////

// C++ BS_REAL array to Python list
%typemap(out) BS_REAL [ANY] {
    int i;
    $result = PyList_New($1_dim0);
    for (i = 0; i < $1_dim0; i++) {
            PyObject *o = PyFloat_FromDouble((double) $1[i]);
            PyList_SetItem($result, i, o);
    }
}

// C++ OpacityFlags to Python dict
%typemap(out) OpacityFlags * {
    PyObject *out_dict = PyDict_New();
    OpacityFlags * in_var = $1;
    PyDict_SetItemString(out_dict, "use_abs_em", PyInt_FromLong(in_var->use_abs_em)); 
    PyDict_SetItemString(out_dict, "use_pair", PyInt_FromLong(in_var->use_pair)); 
    PyDict_SetItemString(out_dict, "use_brem", PyInt_FromLong(in_var->use_brem)); 
    PyDict_SetItemString(out_dict, "use_inelastic_scatt", PyInt_FromLong(in_var->use_inelastic_scatt)); 
    PyDict_SetItemString(out_dict, "use_iso", PyInt_FromLong(in_var->use_iso)); 
    
    $result = out_dict;
}

// C++ OpacityParams to Python dict
%typemap(out) OpacityParams * {
    PyObject *out_dict = PyDict_New();
    OpacityParams * in_var = $1;
    PyDict_SetItemString(out_dict, "use_dU", PyInt_FromLong(in_var->use_dU)); 
    PyDict_SetItemString(out_dict, "use_dm_eff", PyInt_FromLong(in_var->use_dm_eff)); 
    PyDict_SetItemString(out_dict, "use_WM_ab", PyInt_FromLong(in_var->use_WM_ab)); 
    PyDict_SetItemString(out_dict, "use_WM_sc", PyInt_FromLong(in_var->use_WM_sc)); 
    PyDict_SetItemString(out_dict, "use_decay", PyInt_FromLong(in_var->use_decay)); 
    PyDict_SetItemString(out_dict, "use_BRT_brem", PyInt_FromLong(in_var->use_BRT_brem)); 
    PyDict_SetItemString(out_dict, "use_NN_medium_corr", PyInt_FromLong(in_var->use_NN_medium_corr)); 
    PyDict_SetItemString(out_dict, "neglect_blocking", PyInt_FromLong(in_var->neglect_blocking)); 
    
    $result = out_dict;
}

// C++ SpectralOpacities to Python dict
%typemap(out) SpectralOpacities {
    PyObject *j = PyList_New(total_num_species);
    PyObject *j_s = PyList_New(total_num_species);
    PyObject *kappa = PyList_New(total_num_species);
    PyObject *kappa_s = PyList_New(total_num_species);

    PyObject *out_dict = PyDict_New();

    int i;
    SpectralOpacities *in_var = &$1;

    for (i = 0; i < total_num_species; i++) {
            PyList_SetItem(j, i, PyFloat_FromDouble((double) in_var->j[i]));
    }
    PyDict_SetItemString(out_dict, "j", j);

    for (i = 0; i < total_num_species; i++) {
            PyList_SetItem(j_s, i, PyFloat_FromDouble((double) in_var->j_s[i]));
    }
    PyDict_SetItemString(out_dict, "j_s", j_s);

    for (i = 0; i < total_num_species; i++) {
            PyList_SetItem(kappa, i, PyFloat_FromDouble((double) in_var->kappa[i]));
    }
    PyDict_SetItemString(out_dict, "kappa", kappa);

    for (i = 0; i < total_num_species; i++) {
            PyList_SetItem(kappa_s, i, PyFloat_FromDouble((double) in_var->kappa_s[i]));
    }
    PyDict_SetItemString(out_dict, "kappa_s", kappa_s);

    $result = out_dict;
}

// C++ M1Opacities to Python dict
%typemap(out) M1Opacities {
    PyObject *eta = PyList_New(total_num_species);
    PyObject *eta_0 = PyList_New(total_num_species);
    PyObject *kappa_a = PyList_New(total_num_species);
    PyObject *kappa_0_a = PyList_New(total_num_species);
    PyObject *kappa_s = PyList_New(total_num_species);

    PyObject *out_dict = PyDict_New();
    
    int i;
    M1Opacities *in_var = &$1;

    for (i = 0; i < total_num_species; i++) {
            PyList_SetItem(eta, i, PyFloat_FromDouble((double) in_var->eta[i]));
    }
    PyDict_SetItemString(out_dict, "eta", eta);        
    
    for (i = 0; i < total_num_species; i++) {
            PyList_SetItem(eta_0, i, PyFloat_FromDouble((double) in_var->eta_0[i]));
    }
    PyDict_SetItemString(out_dict, "eta_0", eta_0);        

    for (i = 0; i < total_num_species; i++) {
            PyList_SetItem(kappa_a, i, PyFloat_FromDouble((double) in_var->kappa_a[i]));
    }
    PyDict_SetItemString(out_dict, "kappa_a", kappa_a);        

    for (i = 0; i < total_num_species; i++) {
            PyList_SetItem(kappa_0_a, i, PyFloat_FromDouble((double) in_var->kappa_0_a[i]));
    }
    PyDict_SetItemString(out_dict, "kappa_0_a", kappa_0_a);        

    for (i = 0; i < total_num_species; i++) {
            PyList_SetItem(kappa_s, i, PyFloat_FromDouble((double) in_var->kappa_s[i]));
    }
    PyDict_SetItemString(out_dict, "kappa_s", kappa_s);        

    $result = out_dict;
}

// Define Python code for printing output easily
%pythoncode %{
def print_reactions(opacity_flags):
    print('# Included reactions:')
    print('# Neutrino abs on nucleons - e+/e- capture : %d' %opacity_flags["use_abs_em"])
    print('# Elastic scattering on nucleons           : %d' %opacity_flags["use_iso"])
    print('# Nucleon-nucleon bremmstrahlung           : %d' %opacity_flags["use_brem"])
    print('# e+ e- annihilation (and inverse)         : %d' %opacity_flags["use_pair"])
    print('# Inelastic scattering on e+/e-            : %d' %opacity_flags["use_inelastic_scatt"])
    print("")

def print_corrections(opacity_pars):
    print('# Included corrections:')
    print('# dU correction                        : %d' %opacity_pars["use_dU"])
    print('# Effective mass correction            : %d' %opacity_pars["use_dm_eff"])
    print('# Weak magnetism on charged reactions  : %d' %opacity_pars["use_WM_ab"])
    print('# Weak magnetism on neutral reactions  : %d' %opacity_pars["use_WM_sc"])
    print('# Include (inverse) nucleon decays     : %d' %opacity_pars["use_decay"])
    print('# Neglect blocking factors             : %d' %opacity_pars["neglect_blocking"])
    print('# Bremsstrahlung as in Burrows+2006    : %d' %opacity_pars["use_BRT_brem"])
    print('# Medium Fischer+16 correction to brem : %d' %opacity_pars["use_NN_medium_corr"])
    print("")

def print_neutrino_quantities(m1_pars):
    print("    nue            anue            nux            anux");
    print("n  %13.6e  %13.6e   %13.6e  %13.6e     (cm^-3)" %(m1_pars.n[0]  , m1_pars.n[1]  , m1_pars.n[2]  , m1_pars.n[3]))
    print("J  %13.6e  %13.6e   %13.6e  %13.6e (MeV cm^-3)" %(m1_pars.n[0]  , m1_pars.J[1]  , m1_pars.J[2]  , m1_pars.J[3]))
    print("chi %12.10f   %12.10f    %12.10f   %12.10f"   %(m1_pars.chi[0], m1_pars.chi[1], m1_pars.chi[2], m1_pars.chi[3]))
    print("")

def print_spectral_rates(spec_rates):
    print("      j             j_s           kappa         kappa_s");
    print(" nue %13.6e %13.6e %13.6e %13.6e" %(spec_rates['j'][0], spec_rates['j_s'][0], spec_rates['kappa'][0], spec_rates['kappa_s'][0]))
    print("anue %13.6e %13.6e %13.6e %13.6e" %(spec_rates['j'][1], spec_rates['j_s'][1], spec_rates['kappa'][1], spec_rates['kappa_s'][1]))
    print(" nux %13.6e %13.6e %13.6e %13.6e" %(spec_rates['j'][2], spec_rates['j_s'][2], spec_rates['kappa'][2], spec_rates['kappa_s'][2]))
    print("anux %13.6e %13.6e %13.6e %13.6e" %(spec_rates['j'][3], spec_rates['j_s'][3], spec_rates['kappa'][3], spec_rates['kappa_s'][3]))
    print("")

def print_integrated_rates(gray_rates):
    print("      eta0          eta1          kappa0        kappa1        scat1");
    print(" nue %13.6e %13.6e %13.6e %13.6e %13.6e" %(gray_rates['eta_0'][0], gray_rates['eta'][0], gray_rates['kappa_0_a'][0], gray_rates['kappa_a'][0], gray_rates['kappa_s'][0]))
    print("anue %13.6e %13.6e %13.6e %13.6e %13.6e" %(gray_rates['eta_0'][1], gray_rates['eta'][1], gray_rates['kappa_0_a'][1], gray_rates['kappa_a'][1], gray_rates['kappa_s'][1]))
    print(" nux %13.6e %13.6e %13.6e %13.6e %13.6e" %(gray_rates['eta_0'][2], gray_rates['eta'][2], gray_rates['kappa_0_a'][2], gray_rates['kappa_a'][2], gray_rates['kappa_s'][2]))
    print("anux %13.6e %13.6e %13.6e %13.6e %13.6e" %(gray_rates['eta_0'][3], gray_rates['eta'][3], gray_rates['kappa_0_a'][3], gray_rates['kappa_a'][3], gray_rates['kappa_s'][3]))
    print("")
%}

// BNS_NURATES header files
%include "bns_nurates.hpp"
%include "constants.hpp"
%include "distribution.hpp"
%include "eos.hpp"
%include "formulas_beta_iso.hpp"
%include "functions.hpp"
%include "integration.hpp"
%include "kernel_brem.hpp"
%include "kernel_nes.hpp"
%include "kernel_pair.hpp"
%include "kernels.hpp"
%include "m1_opacities.hpp"
%include "weak_magnetism.hpp"
