import numpy as np
import util
from netCDF4 import Dataset as netcdf

class Database:
   def get(self, index):
      values = self._data[:,:,index]
      return values

   # Returns a (weighted) random segment
   #@profile
   def get_random(self, target_state, metric):
      weights = np.zeros(self.size(), float)
      weights = 1.0/metric.compute(target_state, self._data[0,:,:])

      # Do a weighted random choice of the weights
      I = util.random_weighted(weights)
      return self.get(I)


# Trajectories based on gaussian random walk
class Random(Database):
   def __init__(self, N, T, V, variance=1):
      self._N = N
      self._T = T
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
   def __init__(self, filename):
      self._filename = filename
      self._file = netcdf(self._filename)
      vars = self._file.variables
      self._vars = [var for var in vars if var not in ["date", "leadtime"]]
      self._vars = ["air_temperature_2m"]
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
