import numpy as np
# Must implement:
#   Compute how close the two states are
#   def compute(self, state1, state2):
class Metric:
   def __init__(self):
      pass

class Rmsd(Metric):
   @staticmethod
   def compute(state1, state2):
      total = 0
      for var in state1:
         curr = (state1[var] - state2[var])**2
         total = total + curr
      return np.sqrt(total)

class Weighted(Metric):
   # weights: a dictionary of variable:weight entries
   def __init__(self, weights):
      self._weights = weights

   def compute(state1, state2):
      total = 0
      for var in state1:
         weight = self._weights[var]
         curr = weight * (state1[var] - state2[var])**2
         total = total + curr
      return np.sqrt(total)
