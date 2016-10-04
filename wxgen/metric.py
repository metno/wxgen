import numpy as np
# Must implement:
#   Compute how close the two states are
#   def compute(self, state1, state2):
class Metric:
   def __init__(self):
      pass

class Rmsd(Metric):
   @staticmethod
   #@profile
   def compute(state1, state2):
      if state1.shape != state2.shape:
         total = np.sum(abs(np.resize(state1, state2.shape) - state2)**2, axis=0)
      else:
         total = np.sum(state1 - state2)**2
      return np.sqrt(total)

class Weighted(Metric):
   # weights: a dictionary of variable:weight entries
   def __init__(self, weights):
      self._weights = weights

   def compute(state1, state2):
      total = np.sum(weights * (state1 - state2)**2)
      return np.sqrt(total)

class Exp(Metric):
   @staticmethod
   #@profile
   def compute(state1, state2):
      if state1.shape != state2.shape:
         total = np.sum(abs(np.resize(state1, state2.shape) - state2)**2, axis=0)
      else:
         total = np.sum(state1 - state2)**2
      return np.sqrt(total)

