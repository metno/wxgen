import inspect
import numpy as np
import sys
import wxgen.util

# Returns a list of all metric classes
def get_all():
   temp = inspect.getmembers(sys.modules[__name__], inspect.isclass)
   return temp


# Returns a metric object of a class with the given name
def get(name):
   metrics = get_all()
   m = None
   for mm in metrics:
      if(name == mm[0].lower()):
         m = mm[1]()
   if m is None:
      wxgen.util.error("Cannot find metric called '%s'" % name)
   return m



# Must implement:
#   Compute how close the two states are
#   def compute(self, state1, state2):
class Metric(object):
   # +1 means higher values are better, -1 means lower values are better
   _orientation = 1

   def __init__(self):
      pass

   def compute(self, state1, state2):
      # Ensure states are the same size
      if state1.shape != state2.shape:
         state1 = wxgen.util.resize(state1, state2.shape)
      return self._compute(state1, state2)


# The score is diff ** 2
class Rmsd(Metric):
   _orientation = -1
   def _compute(self, state1, state2):
      assert(state1.shape == state2.shape)
      total = np.sum(abs(state1 - state2)**2, axis=0)
      return np.sqrt(total)

# THe score is (weights * diff)**2
class Weighted(Metric):
   _orientation = -1
   # weights: an array of variable-weights
   def __init__(self, weights):
      self._weights = weights

   def _compute(self, state1, state2):
      if self._weights is None:
         weights = 1
      else:
         weights = wxgen.util.resize(self._weights, state2.shape)
      total = np.sum(weights*abs(state1 - state2)**2, axis=0)
      return np.sqrt(total)

# The score is exp(-factor * diff)
class Exp(Metric):
   _orientation = 1
   def __init__(self, factors=None):
      self._factors = factors

   def _compute(self, state1, state2):
      if self._factors is None:
         factors = 1
      else:
         factors = wxgen.util.resize(self._factors, state2.shape)
      total = np.exp(np.sum(-factors*abs(state1 - state2), axis=0))
      return total

