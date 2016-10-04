import numpy as np
import util
from netCDF4 import Dataset as netcdf

class Database(object):
   def get(self, index):
      values = self._data[:,:,index]
      return values

   # Returns a (weighted) random segment
   #@profile
   def get_random(self, target_state, metric):
      weights = metric.compute(target_state, self._data[0,:,:])

      # Flip the metric if it is negative oriented
      if metric._orientation == -1:
         I0 = np.where(weights < 1e-3)[0]
         I1 = np.where(weights >= 1e-3)[0]
         # Ensure we do not get too high weights
         weights[I1] = 1.0/weights[I1]
         weights[I0] = 1e3

      # Do a weighted random choice of the weights
      I = util.random_weighted(weights)
      return self.get(I)


# Trajectories based on gaussian random walk
class Random(Database):
   def __init__(self, N, T, V, variance=1):
      self._N = N
      self._T = T
      if V == None:
         V = 1
      self._V = V
      self._variance = variance
      self._data = np.zeros([T, V, N], float)
      for v in range(0, self._V):
         self._data[:,v,:]  = np.cumsum(np.random.randn(T, N)*np.sqrt(self._variance), axis=0)

   # Number of days
   def days(self):
      return self._T

   # Number of trajectories
   def size(self):
      return self._N

   def vars(self):
      return range(0, self._V)

   def num_vars(self):
      return self._V


class Netcdf(Database):
   def __init__(self, filename, V=None):
      self._filename = filename
      self._file = netcdf(self._filename)
      vars = self._file.variables
      self._vars = [var for var in vars if var not in ["date", "leadtime"]]
      if V is not None:
         self._vars = self._vars[0:V]
      self._size = self._num_members() * self._file.dimensions["date"].size
      self._num_vars = len(self._vars)

      # Load data
      V = self._num_vars
      N = self.size()
      T = self.days()
      self._data = np.zeros([T, V, N], float)
      for v in range(0, V):
         var = self.vars()[v]
         temp = self._file.variables[var]
         for t in range(0, T):
            self._data[t,v,:] = temp[:,t,:].flatten()

   # Number of days
   def days(self):
      return self._file.dimensions["leadtime"].size

   # Number of trajectories
   def size(self):
      return self._size

   def _num_members(self):
      return self._file.dimensions["member"].size

   def vars(self):
      return self._vars

   def num_vars(self):
      return self._num_vars
