import numpy as np
import wxgen.util
import datetime
import wxgen.variable
import copy
try:
   from netCDF4 import Dataset as netcdf
except:
   from scipy.io.netcdf import netcdf_file as netcdf


class Database(object):
   """
   Abstract class which stores weather segments (short time-series)

   The database is stored in self._data, which has dimensions (T, V, N, X, Y) where T is the segment
   length, V is the number of variables, N is the number of segments, and X and Y are geographical
   dimensions. A subclass must populate this field during initialization.

   Attributes:
   variables      A list of wxgen.variable.Variable

   """
   def __init__(self):
      self._debug = False
      """
      External state (such as the month of the year)
      """
      self._ext_state = None

   def info(self):
      print "Database information:"
      print "  Length of segments: %d" % self.days()
      print "  Number of segments: %d" % self.size()
      print "  Number of variables: %d" % len(self.variables)

   def size(self):
      """ Get the number of segments in the database """
      raise NotImplementedError()

   def days(self):
      """ Get the length of a segment in the database """
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

      self.variables = [wxgen.variable.Variable(str(i)) for i in range(0, self._V)]

   def days(self):
      return self._T

   def size(self):
      return self._N


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
      self.variables = [wxgen.variable.Variable(name) for name in self._file.variables if name not in ["date", "leadtime"]]
      if V is not None:
         self.variables = self.variables[0:V]
      self._size = self._num_members() * self._file.dimensions["date"].size

      # Load data
      V = len(self.variables)
      N = self.size()
      T = self.days()
      self._data = np.zeros([T, V, N], float)
      self._ext_state = np.zeros(N, float)
      assert(self._ext_state.shape[0] == N)
      for v in range(0, V):
         var = self.variables[v]
         temp = self._copy(self._file.variables[var.name]) # dims: date,leadtime,member

         # Quality control
         if var.name == "precipitation_amount":
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

   def _copy(self, data):
      data = data[:].astype(float)
      q = copy.deepcopy(data)
      # Remove missing values. Convert to -999 and then back to nan to avoid
      # warning messages when doing <, >, and == comparisons with nan.
      q[np.isnan(q)] = -999
      q[(q == -999) | (q < -1000000) | (q > 1e30)] = np.nan
      return q

class Lorenz63(Database):
   """
   Trajectories based on the Lorenz 63 model.
   """
   def __init__(self, N, T, R=28, S=10, B=2.6667, dt=0.0001):
      Database.__init__(self)
      self._N = N
      self._T = T
      self._R = R
      self._S = S
      self._B = B
      self._dt = dt
      self._V = 3  # Number of variables
      self._data = np.zeros([T, self._V, N], float)
      self._initial_state = [-10, -10, 25]
      self._std_initial_state = 0.1  # Standard deviation of initial condition error

      # Initialize
      for v in range(0, self._V):
         self._data[0,v,:] = self._initial_state[v] + np.random.randn(N) * self._std_initial_state

      TT = int(1 / self._dt)/10
      # Iterate
      for t in range(1, T):
         x0 = copy.deepcopy(self._data[t-1,0,:])
         y0 = copy.deepcopy(self._data[t-1,1,:])
         z0 = copy.deepcopy(self._data[t-1,2,:])
         # Iterate from t-1 to t
         for tt in range(0, TT):
            x1 = x0 + self._dt * self._S * (y0 - x0)
            y1 = y0 + self._dt * (self._R * x0 - y0 - (x0 * z0))
            z1 = z0 + self._dt * (x0*y0 - self._B*z0)

            x2 = x1 + self._dt * self._S * (y1-x1)
            y2 = y1 + self._dt * (self._R * x1 - y1 - x1*z1)
            z2 = z1 + self._dt * (x1*y1 - self._B*z1)

            x0 = 0.5 * (x2 + x0)
            y0 = 0.5 * (y2 + y0)
            z0 = 0.5 * (z2 + z0)

         self._data[t,0,:] = x0
         self._data[t,1,:] = y0
         self._data[t,2,:] = z0

      self.variables = [wxgen.variable.Variable(i) for i in ["X", "Y", "Z"]]

   def days(self):
      return self._T

   def size(self):
      return self._N
