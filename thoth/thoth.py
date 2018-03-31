import numbers
from _bstrap_bc import * 
import numpy as np

def calc_entropy(array, n_samp):
   types = (list, np.ndarray, tuple)
   if (not isinstance(array, types)):
      print("Error! Please input an array")
      print("Support types are: list, numpy array, and tuple")
      return
   else:
      if all(isinstance(x, numbers.Number) for x in array):
         prob = prob_from_array(array)
         norm_prob(prob)
         output = entropy_bs_wrapper(prob, n_samp, rng(), 5)
      else:
         print("Error! Not all values in array are numeric")
         return
   return output

def calc_jsd(array1, array2, alpha, n_samp):
   types = (list, np.ndarray, tuple)
   if ((not isinstance(array1, types)) or (not isinstance(array2, types))):
      print("Error! Please input two arrays")
      print("Support types are: list, numpy array, and tuple")
      return
   else:
      if (all(isinstance(x, numbers.Number) for x in array1) and (all(isinstance(x, numbers.Number) for x in array2))):
         prob1 = prob_from_array(array1)
         prob2 = prob_from_array(array2)
         norm_prob(prob1)
         norm_prob(prob2)
         output = jsd_bs_wrapper(prob1, prob2, alpha, n_samp, rng(), 5)
      else:
         print("Error! Not all values in the arrays are numeric")
         return
   return output

def calc_mi(arrays, n_samp):
   types = (list, np.ndarray, tuple)
   check_arrays = isinstance(arrays, types) and all(isinstance(x, types) for x in arrays)
   if not check_arrays:
      print("Error! Please an array of arrays")
      print("Support types are: list, numpy array, and tuple")
      return
   else:
      if all(isinstance(x, numbers.Number) for y in arrays for x in y):
         j_prob = jprob_from_array(arrays)
         norm_jprob(j_prob)
         output = mi_bs_wrapper(j_prob, n_samp, rng(), 5)
      else:
         print("Error! Not all values in the arrays are numeric")
         return
   return output
