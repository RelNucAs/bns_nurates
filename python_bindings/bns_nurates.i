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
%typemap(in) struct OpacityParams *OpacityParams::OpacityParams {
	$1.use_dU = PyDict_GetItemString($input, 'use_dU');
	$1.use_dm_eff = PyDict_GetItemString($input, 'use_dm_eff');
	$1.use_WM_ab = PyDict_GetItemString($input, 'use_WM_ab');
	$1.use_WM_sc = PyDict_GetItemString($input, 'use_WM_sc');
	$1.use_decay = PyDict_GetItemString($input, 'use_decay');
	$1.neglect_blocking = PyDict_GetItemString($input, 'neglect_blocking');
}
	

%typemap(in) struct OpacityFlags * {
	PyObject *o = PyDict_GetItemString($input, "use_abs_em");
	printf("%d\n", (int) PyInt_AsLong(o));
	$1->use_abs_em = (int) PyLong_AsLong(o);
	o = PyDict_GetItemString($input, "use_pair");
	$1->use_pair = (int) PyLong_AsLong(o);
	o = PyDict_GetItemString($input, "use_brem");
	$1->use_brem = (int) PyLong_AsLong(o);
	o = PyDict_GetItemString($input, "use_inelastic_scatt");
	$1->use_inelastic_scatt = (int) PyLong_AsLong(o);
	o = PyDict_GetItemString($input, "use_iso");
	$1->use_iso = (int) PyLong_AsLong(o);
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
