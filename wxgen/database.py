import numpy as np
import util
try:
   from netCDF4 import Dataset as netcdf
except:
   from scipy.io.netcdf import netcdf_file as netcdf


class Database(object):
   """
   Abstract class which stores weather segments (short time-series)
   """
   def __init__(self):
      self._num_vars = None
      self._debug = False

   def info(self):
      print "Database information:"
      print "  Length of segments: %d" % self.days()
      print "  Number of segments: %d" % self.size()
      print "  Number of variables: %d" % self.num_vars()

   def num_vars(self):
      """ Get the number of variables in the database """
      if self._num_vars is None:
         self._num_vars = len(self.vars())
      return self._num_vars

   def size(self):
      """ Get the number of segments in the database """
      raise NotImplementedError()

   def days(self):
      """ Get the length of a segment in the database """
      raise NotImplementedError()

   def vars(self):
      """ Get the names of the variables in the database """
      raise NotImplementedError()

   def get(self, i):
      """ Get the i'th trajectory in the database """
      values = self._data[:,:,i]
      return values

   def get_random(self, target_state, metric):
      """
      Returns a random segment from the database that is weighted
      by the scores computed by metric.

      target_state   A numpy array (length V)
      metric         Of type wxgen.metric.Metric
      """
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
      if self._debug:
         print "Data: ", self._data[0,:,I]
         print "Weight: ", weights[I]
         print "Max weight: ", np.max(weights)
      return self.get(I)


class Random(Database):
   """
   Trajectories are random Gaussian walks, with constant variance over time. There is no covariance
   between different forecast variables.
   """
   def __init__(self, N, T, V, variance=1):
      Database.__init__()
      self._N = N
      self._T = T
      if V == None:
         V = 1
      self._V = V
      self._variance = variance
      self._data = np.zeros([T, V, N], float)

      # Ensure that the signal has a constant variance over time
      scale = 1./np.sqrt(np.linspace(1, T, T))

      for v in range(0, self._V):
         self._data[:,v,:]  = np.transpose(np.resize(scale, [N, T])) * np.cumsum(np.random.randn(T, N)*np.sqrt(self._variance), axis=0)

   def days(self):
      return self._T

   def size(self):
      return self._N

   def vars(self):
      return range(0, self._V)


class Netcdf(Database):
   """
   Segments stored in a netcdf database

   Should have the following format:
      dims: date, leadtime, member
      vars: date(date), leadtime(leadtime)
            variable_name(date, leadtime, member)
   where variable_name is one or more names of weather variables
   """
   def __init__(self, filename, V=None):
      """
      filename    Load data from this file
      V           Only use the first V variables in the database
      """
      Database.__init__(self)
      self._file = netcdf(filename)

      # Set dimensions
      self._vars = [var for var in self._file.variables if var not in ["date", "leadtime"]]
      if V is not None:
         self._vars = self._vars[0:V]
      self._size = self._num_members() * self._file.dimensions["date"].size

      # Load data
      V = self.num_vars()
      N = self.size()
      T = self.days()
      self._data = np.zeros([T, V, N], float)
      for v in range(0, V):
         var = self.vars()[v]
         temp = self._file.variables[var]
         for t in range(0, T):
            self._data[t,v,:] = temp[:,t,:].flatten()

   def days(self):
      return self._file.dimensions["leadtime"].size

   def size(self):
      return self._size

   def _num_members(self):
      return self._file.dimensions["member"].size

   def vars(self):
      return self._vars
