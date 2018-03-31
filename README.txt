===================
*      THOTH      *
===================
****version 1.0****
===================

/* Once you've installed: see the test file 'test_v1.py' for some demo usage */

/* Installation Guide */

REQUIREMENTS:
Python-2.7.3
GSL-1.15
SWIG-2.0 (If you want to write your own functions)

0. This installation is known to work easily if you use Mac OS X and have 
installed python and gsl using MacPorts.

1. Modify the setup_cfg.py
Change the 'GSL_INCLUDE' and 'GSL_LIB' variables to your path to the GSL
header and library files. For example:
GSL_INCLUDE = 'some_path/gsl-1.15/include'
GSL_LIB = 'some_path/gsl-1.15/lib'


2. Then open up a terminal and type
user_name$ python setup.py install -dir-

replace '-dir-' for the directory you want the package to be installed.
For example: user_name$ python setup.py install /home/user_name/thoth

Without specifying '-dir-', the package will be installed into the 
default python directory so that you can simply import in in Python
without adding to 'sys.path'.


/* Load the package in Python */
If you have GSL installed in a custom directory, you need to modify your LD_LIBRARY_PATH to
include the GSL shared libraries. For example, in a bash shell:
LD_LIBRARY_PATH=some_path/gsl-1.15/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

Then in Python:
If you install the package into Python default directory:
>>> import thoth.thoth as thoth
>>> dir(thoth)

If you install the package into a custom directory, say /home/user_name/thoth/ for example:
>>> import sys
>>> sys.path.append('/home/user_name/thoth/')
>>> import thoth.thoth as thoth
>>> dir(thoth)


/* Available functions */

1. prob_from_array(array)
   This function return a python object which store a pointer to the C-struct 'prob' from a python array.

2. jprob_from_array(arrays)
   Similar to 'prob_from_array', takes an array of arrays and return a pointer to 'j_prob'

3. norm_prob(prob *a)
   Normalize prob a.

4. norm_jprob(j_prob *a)
   Normalize j_prob a.

5. print_prob(prob *a)
   Print prob a to screen.

6. print_jprob(j_prob *a)
   Print j_prob a to screen.

7. entropy(prob* a)
   Calculate entropy of prob a.

8. entropy_ww(prob* a)

9. mi(j_prob *a)
   Calculate mutual information of j_prob a.

10. mi_ww(j_prob *pij, double beta)

11. jsd(prob *a, prob *b, double alpha)

12. entropy_nsb(prob *a)
    Calculate NSB for entropy of prob a.

13. mi_nsb(j_prob *a)
    Calculate NSB for mutual information of j_prob a.

14. calc_entropy(array, n_sample)
    Calculate bootstrap for entropy of a probability.

15. calc_jsd(array1, array2, alpha, n_sample)
    Calculate bootstrap for JSD of a probability.

16. calc_mi(arrays, n_sample)
    Calculate bootstrap for mutual information of a probability.

/* Contributing to THOTH */

You can write your own C codes and add it to THOTH. If you don't want to
expose your functions to the user in python, you don't have to care about 
this section.

If you have, for example, a function 'new_add_function'

int new_add_function(int a, int b) {
   return a+b;
}

that you want the user to use in Python, you should add this function to the header file
'prob_api.h'. This file contains the headers of all the C functions available to the user 
in Python.

Then you use SWIG to generate the Python module from the SWIG interface file located
in the 'src' folder:
swig -python bstrap_bc.i

This will replace the C file called 'bstrap_bc_wrap.c' and include your new function.
If you reinstall the package, you should be able to use your new function in python.

Note: for functions that use and return generic data types (no pointers), you don't need
to modify the SWIG interface file 'bstrap_bc.i'. If you want to use pointers, arrays, etc,
please see SWIG documentation for more details.



