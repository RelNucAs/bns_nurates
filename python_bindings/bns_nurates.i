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

%include "../include/bns_nurates.h"
%include "../include/constants.h"
%include "../include/distribution.h"
%include "../include/functions.h"
%include "../include/integration.h"
%include "../include/kernels.h"
%include "../include/opacities.h"
%include "../include/weak_magnetism.h"

%include "numpy.i"

%init %{
import_array();
%}

%include "carrays.i"
%array_class(int, intArray);

