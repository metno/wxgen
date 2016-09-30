import numpy as np
class Database:
   pass

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

   # Number of days
   def days(self):
      pass

   # Number of trajectories
   def size(self):
      pass

   def vars(self):
      pass

   def get(self, index):
      pass

