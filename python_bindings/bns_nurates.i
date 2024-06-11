%module bns_nurates
%{
#define SWIG_FILE_WITH_INIT
#include "../include/bns_nurates.h"
#include "../include/constants.h"
#include "../include/distribution.h"
#include "../include/functions.h"
#include "../include/integration.h"
#include "../include/kernels.h"
#include "../include/opacities.h"
#include "../include/weak_magnetism.h"
%}

#define __attribute__(x) 
%typemap(out) double j[ANY], double kappa[ANY], double j_s[ANY], double kappa_s[ANY] {
        int i;
        $result = PyList_New($1_dim0);
        for (i = 0; i < $1_dim0; i++) {
                PyObject *o = PyFloat_FromDouble((double) $1[i]);
                PyList_SetItem($result, i, o);
        }
}

%typemap(out) double eta_0[ANY], double eta[ANY], double kappa_0_a[ANY], double kappa_a[ANY], double kappa_s[ANY] {
        int i;
        $result = PyList_New($1_dim0);
        for (i = 0; i < $1_dim0; i++) {
                PyObject *o = PyFloat_FromDouble((double) $1[i]);
                PyList_SetItem($result, i, o);
        }
}

%typemap(out) double integrand[ANY] {
        int i;
        $result = PyList_New($1_dim0);
        for (i = 0; i < $1_dim0; i++) {
                PyObject *o = PyFloat_FromDouble((double) $1[i]);
                PyList_SetItem($result, i, o);
        }
}

%typemap(out) double n[ANY], double J[ANY], double chi[ANY] {
        int i;
        $result = PyList_New($1_dim0);
        for (i = 0; i < $1_dim0; i++) {
                PyObject *o = PyFloat_FromDouble((double) $1[i]);
                PyList_SetItem($result, i, o);
        }
}

%typemap(out) double points[ANY], double w[ANY] {
        int i;
        $result = PyList_New($1_dim0);
        for (i = 0; i < $1_dim0; i++) {
                PyObject *o = PyFloat_FromDouble((double) $1[i]);
                PyList_SetItem($result, i, o);
        }
}



%typemap(in) double j[ANY] (double temp[$1_dim0]) {
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
                        temp[i] = (double) PyFloat_AsDouble(o);
                } else {
                        PyErr_SetString(PyExc_ValueError, "Sequence elements must be numbers");
                        SWIG_fail;
                }
        }
        $1 = temp;
}

%typemap(in) double n[ANY] (double temp[$1_dim0]) {
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
                        temp[i] = (double) PyFloat_AsDouble(o);
                } else {
                        PyErr_SetString(PyExc_ValueError, "Sequence elements must be numbers");
                        SWIG_fail;
                }
        }
        $1 = temp;
}

%typemap(in) double J[ANY] (double temp[$1_dim0]) {
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
                        temp[i] = (double) PyFloat_AsDouble(o);
                } else {
                        PyErr_SetString(PyExc_ValueError, "Sequence elements must be numbers");
                        SWIG_fail;
                }
        }
        $1 = temp;
}

%typemap(in) double chi[ANY] (double temp[$1_dim0]) {
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
                        temp[i] = (double) PyFloat_AsDouble(o);
                } else {
                        PyErr_SetString(PyExc_ValueError, "Sequence elements must be numbers");
                        SWIG_fail;
                }
        }
        $1 = temp;
}
/*
%typemap(in) OpacityParams * (OpacityParams tmp) {
        tmp.use_dU = PyDict_GetItemString($input, "use_dU");
	tmp.use_dm_eff = PyDict_GetItemString($input, "use_dm_eff");
	tmp.use_WM_ab = PyDict_GetItemString($input, "use_WM_ab");
	tmp.use_WM_sc = PyDict_GetItemString($input, "use_WM_sc");
	tmp.use_decay = PyDict_GetItemString($input, "use_decay");
	tmp.use_BRT_brem = PyDict_GetItemString($input, "use_BRT_brem");
        tmp.use_NN_medium_corr = PyDict_GetItemString($input, "use_NN_medium_corr");
	tmp.neglect_blocking = PyDict_GetItemString($input, "neglect_blocking");

	$1 = &tmp;
}
*/

%typemap(in) OpacityFlags * (OpacityFlags tmp) {
	tmp.use_abs_em = (int) PyInt_AsLong(PyDict_GetItemString($input, "use_abs_em"));
	tmp.use_pair = (int) PyInt_AsLong(PyDict_GetItemString($input, "use_pair"));
	tmp.use_brem = (int) PyInt_AsLong(PyDict_GetItemString($input, "use_brem"));
	tmp.use_inelastic_scatt = (int) PyInt_AsLong(PyDict_GetItemString($input, "use_inelastic_scatt"));
	tmp.use_iso = (int) PyInt_AsLong(PyDict_GetItemString($input, "use_iso"));
	
	$1 = &tmp;
}

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

%typemap(out) SpectralOpacities ComputeSpectralOpacitiesNotStimulatedAbs, SpectralOpacities ComputeSpectralOpacitiesStimulatedAbs {
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

%typemap(out) M1Opacities ComputeM1Opacities, M1Opacities ComputeM1OpacitiesNotStimulated {
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

%pythoncode %{
def print_reactions(opacity_flags):
	print('# List of included reactions:')
	print('# Neutrino abs on nucleons - e+/e- capture : %d' %opacity_flags["use_abs_em"])
	print('# Elastic scattering on nucleons           : %d' %opacity_flags["use_iso"])
	print('# Nucleon-nucleon bremmstrhalung           : %d' %opacity_flags["use_brem"])
	print('# e+ e- annihilation (and inverse)         : %d' %opacity_flags["use_pair"])
	print('# Inelastic scattering on e+/e-            : %d' %opacity_flags["use_inelastic_scatt"])
	print("")

def print_corrections(opacity_pars):
	print('# List of included corrections:')
	print('# dU correction                       : %d' %opacity_pars["use_dU"])
	print('# Effective mass correction           : %d' %opacity_pars["use_dm_eff"])
	print('# Weak magnetism on charged reactions : %d' %opacity_pars["use_WM_ab"])
	print('# Weak magnetism on neutral reactions : %d' %opacity_pars["use_WM_sc"])
	print('# Neglect blocking                    : %d' %opacity_pars["neglect_blocking"])
	print("")

def print_spectral_rates(spec_rates):
	print("Spectral rates for nue, anue, nux, anux neutrinos:")
	print("j   [s-1]      : ", *spec_rates["j"])
	print("j_s [s-1]      : ", *spec_rates["j_s"])
	print("kappa   [cm-1] : ", *spec_rates["kappa"])
	print("kappa_s [cm-1] : ", *spec_rates["kappa_s"])
	print("")

def print_integrated_rates(int_rates):
	print("Integrated rates for nue, anue, nux, anux neutrinos:")
	print("eta   [cm-3 s-1]     : ", *int_rates["eta"])
	print("eta_0 [MeV cm-3 s-1] : ", *int_rates["eta_0"])
	print("kappa_a   [cm-1]     : ", *int_rates["kappa_a"])
	print("kappa_0_a [cm-1]     : ", *int_rates["kappa_0_a"])
	print("kappa_s   [cm-1]     : ", *int_rates["kappa_s"])
	print("")
%}

/*

%typemap(in) struct GreyOpacityParams * (GreyOpacityPApacityParams tmp) {
        tmp.use_dm_eff = PyDict_GetItemString($input, "use_dm_eff");
        tmp.use_WM_ab = PyDict_GetItemString($input, "use_WM_ab");
        tmp.use_WM_sc = PyDict_GetItemString($input, "use_WM_sc");
        tmp.use_decay = PyDict_GetItemString($input, "use_decay");
        tmp.neglect_blocking = PyDict_GetItemString($input, "neglect_blocking");
        
        $1 = &tmp;
}
*/
/*
%typemap(in) struct OpacityParams *OpacityParams::OpacityParams {
	$1.use_dU = 1.; //PyDict_GetItemString($input, 'use_dU');
	$1.use_dm_eff = PyDict_GetItemString($input, 'use_dm_eff');
	$1.use_WM_ab = PyDict_GetItemString($input, 'use_WM_ab');
	$1.use_WM_sc = PyDict_GetItemString($input, 'use_WM_sc');
	$1.use_decay = PyDict_GetItemString($input, 'use_decay');
	$1.neglect_blocking = PyDict_GetItemString($input, 'neglect_blocking');
}
*/	


//%apply double{j[ANY]};
//%apply (double SpectralOpacities::j[ANY]) {(double SpectralOpacities::j[ANY])};

//%include "carrays.i"
//%array_class(int, intArray);

%include "../include/bns_nurates.h"
%include "../include/constants.h"
%include "../include/distribution.h"
%include "../include/functions.h"
%include "../include/integration.h"
%include "../include/kernels.h"
%include "../include/opacities.h"
%include "../include/weak_magnetism.h"

//%include "numpy.i"

//%init %{
//import_array();
//%}

%pythoncode %{
def PySpectralOpacities(e_nu, quad_1d, grey_pars):
	return ComputeSpectralOpacitiesNotStimulatedAbs(e_nu, quad_1d, grey_pars)
	
%}
