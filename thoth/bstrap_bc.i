%module bstrap_bc
%{
#define SWIG_FILE_WITH_INIT
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include "prob.h"
%}

%include "numpy.i"

%init %{
	import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(double* in_array, int n)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* j_in_array, int n1, int n2)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* out_results, int quantities)}

%rename(create_new_prob) new_prob(int n);
%rename(remove_delete_prob) delete_prob(prob *a);

%include "prob_api.h"
