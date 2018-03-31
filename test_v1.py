import sys
sys.path.append('/home/quanbui/THOTH/')
import thoth.thoth as thoth
import numpy as np

print("=========ENTROPY===========")
a_prob = thoth.prob_from_array([2,3,12,5,7])
thoth.norm_prob(a_prob)
thoth.entropy(a_prob)
thoth.entropy_nsb(a_prob)
results = thoth.calc_entropy([2,3,12,5,7], 100000)
print(results)

print("=========JSD===========")
results = thoth.calc_jsd([1,2,3,4,5],[5,4,3,2,1], 0.5, 100000)
print(results)

print("=========MI===========")
arrays = np.array([[1,2,3], [3,5,2], [3,2,4]])
results = thoth.calc_mi(arrays, 100000)
print(results)

