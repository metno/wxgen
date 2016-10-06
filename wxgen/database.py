import numpy as np
import wxgen.util
import datetime
import copy
try:
   from netCDF4 import Dataset as netcdf
except:
   from scipy.io.netcdf import netcdf_file as netcdf


class Database(object):
   """
   Abstract class which stores weather segments (short time-series)

   The database is stored in self._data, which has dimensions (T, V, N) where T is the segment
   length, V is the number of variables, and N is the number of segments. A subclass must populate
   this field during initialization.
   """
   def __init__(self):
      self._num_vars = None
      self._debug = False
      """
      External state (such as the month of the year)
      """
      self._ext_state = None

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

   def get_truth(self):
      trajectory = np.squeeze(self._data[0, :, :])
      trajectory = trajectory.transpose()
      return trajectory

   def get_random(self, target_state, metric, ext_state=None):
      """
      Returns a random segment from the database that is weighted
      by the scores computed by metric.

      target_state   A numpy array (length V)
      metric         Of type wxgen.metric.Metric
      ext_state      External state
      """
      weights = metric.compute(target_state, self._data[0,:,:])

      # Flip the metric if it is negative oriented
      if metric._orientation == -1:
         I0 = np.where(weights < 1e-3)[0]
         I1 = np.where(weights >= 1e-3)[0]
         # Ensure we do not get too high weights
         weights[I1] = 1.0/weights[I1]
         weights[I0] = 1e3

      # ext state
      Iall = range(0, len(weights))
      if ext_state is not None and self._ext_state is not None:
         II = np.where(self._ext_state == ext_state)[0]
         if len(II) == 0:
            wxgen.util.error("Cannot find a segment with  external state = %s" % str(ext_state))
         weights = weights[II]
         I = wxgen.util.random_weighted(weights)
         I = II[I]
      else:
         I = wxgen.util.random_weighted(weights)

      # Do a weighted random choice of the weights
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
      Database.__init__(self)
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
   _debug0 = False
   def __init__(self, filename, V=None):
      """
      filename    Load data from this file
      V           Only use the first V variables in the database
      """
      Database.__init__(self)
      self._file = netcdf(filename)

      self._datename = "date"

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
      self._ext_state = np.zeros(N, float)
      assert(self._ext_state.shape[0] == N)
      for v in range(0, V):
         var = self.vars()[v]
         temp = self._copy(self._file.variables[var]) # dims: date,leadtime,member

         # Quality control
         if var == "precipitation_amount":
            temp[temp < 0] = np.nan
         for t in range(0, T):
            self._data[t,v,:] = temp[:,t,:].flatten()
            if self._debug0 and t == 0:
               print self._data[t, v, :]
               print temp[0:2,t,0:2]
               print temp[0:2,t,0]
               print temp[0:2,t,0:2].flatten()
      times = self._file.variables[self._datename]
      day_of_year = np.zeros(times.shape)
      for d in range(0, times.shape[0]):
         day_of_year[d] = datetime.datetime.fromtimestamp(times[d]).strftime('%j')
      month_of_year = day_of_year.astype(int)/ 30
      self._ext_state = np.repeat(month_of_year, self._num_members())
      if self._debug0:
         print self._ext_state


   def days(self):
      return self._file.dimensions["leadtime"].size

   def size(self):
      return self._size

   def _num_members(self):
      return self._file.dimensions["member"].size

   def vars(self):
      return self._vars

   def _copy(self, data):
      data = data[:].astype(float)
      q = copy.deepcopy(data)
      # Remove missing values. Convert to -999 and then back to nan to avoid
      # warning messages when doing <, >, and == comparisons with nan.
      q[np.isnan(q)] = -999
      q[(q == -999) | (q < -1000000) | (q > 1e30)] = np.nan
      return q
