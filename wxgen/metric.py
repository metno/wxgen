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

# The score is diff ** 2
class Rmsd(Metric):
   _orientation = -1
   #@profile
   def compute(self, state1, state2):
      if state1.shape != state2.shape:
         total = np.sum(abs(np.resize(state1, state2.shape) - state2)**2, axis=0)
      else:
         total = np.sum(state1 - state2)**2
      return np.sqrt(total)

# THe score is (weights * diff)**2
class Weighted(Metric):
   _orientation = -1
   # weights: an array of variable-weights
   def __init__(self, weights):
      self._weights = weights

   def compute(self, state1, state2):
      if state1.shape != state2.shape:
         total = np.sum(np.resize(self._weights, state2.shape)*abs(np.resize(state1, state2.shape) - state2)**2, axis=0)
      else:
         total = np.sum(self._weights*(state1 - state2)**2)
      return np.sqrt(total)

# The score is exp(-factor * diff)
class Exp(Metric):
   _orientation = 1
   def __init__(self, factors=None):
      self._factors = factors

   #@profile
   def compute(self, state1, state2):
      F = np.resize(self._factors, state2.shape)
      if state1.shape != state2.shape:
         s1 = np.reshape(np.repeat(state1, state2.shape[1]), state2.shape)
         total = np.exp(np.sum(-F*abs(s1 - state2), axis=0))
      else:
         total = np.exp(np.sum(-F*abs(state1 - state2)))
      return total

