import numpy as np
import util
from netCDF4 import Dataset as netcdf

class Database:
   # Returns the trajectory with the best metric
   def get_min(self, target_state, metric):
      weights = np.zeros(self.size(), float)
      for i in range(0, self.size()):
         curr = self.get(i)
         state = dict()
         for var in self.vars():
            state[var] = curr[var][0]
         weights[i] = 1.0/metric.compute(target_state, state)

      # Do a weighted random choice of the weights
      I = util.random_weighted(weights)
      return self.get(I)


# Trajectories based on gaussian random walk
class Random(Database):
   _variance = 1
   # Number of days
   def days(self):
      return 10

   # Number of trajectories
   def size(self):
      return 3000

   def vars(self):
      return ["T"]

   def get(self, index):
      data = dict()

      T = self.days()
      for var in self.vars():
         #data[var] = np.cumsum(np.random.randn(T))/np.sqrt(range(1, T+1))
         #data[var] = np.cumsum(np.random.randn(T))*np.exp(-0.01*np.linspace(0, T, T))
         #data[var] = np.random.randn(1)*1+ np.cumsum(np.random.randn(T))
         data[var] = np.cumsum(np.random.randn(T)*np.sqrt(self._variance))
      return data


class Netcdf(Database):
   def __init__(self, filename):
      self._filename = filename
      self._file = netcdf(self._filename)
      self._data = dict()
      vars = self._file.variables
      self._vars = [var for var in vars if var not in ["date", "leadtime"]]
      for var in self.vars():
         self._data[var] = self._file.variables[var]

   # Number of days
   def days(self):
      return self._file.dimensions["leadtime"].size

   # Number of trajectories
   def size(self):
      return self._num_members() * self._file.dimensions["date"].size

   def _num_members(self):
      return self._file.dimensions["member"].size

   def vars(self):
      return self._vars

   #@profile
   def get(self, index):
      d = index / self._num_members()
      m = index % self._num_members()
      values = dict()
      for var in self.vars():
         values[var] = self._data[var][d, :, m]
      return values

