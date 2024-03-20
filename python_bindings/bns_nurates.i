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

%typemap(out) double eta_0[ANY], double eta[ANY], double kappa_a_0[ANY], double kappa_a[ANY], double kappa_s[ANY] {
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


