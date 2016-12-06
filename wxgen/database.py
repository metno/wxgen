import numpy as np
import wxgen.util
import datetime
import wxgen.aggregator
import wxgen.variable
import copy
import netCDF4


class Database(object):
   """
   Abstract class which stores weather segments (short time-series)

   The database is stored in self._data, which has dimensions (T, X, Y, V, N) where T is the segment
   length, V is the number of variables, N is the number of segments, and X and Y are geographical
   dimensions. A subclass must populate this field during initialization.

   Attributes:
   variables      A list of wxgen.variable.Variable
   length         Length of each trajectory (in number of days)
   num            Number of trajectories
   lats           2D grid of latitudes
   lons           2D grid of longitudes
   X              Number of X-axis points
   Y              Number of Y-axis points
   inittimes      numpy array with initialization time corresponding to each member

   Internal:
   _data          A 5D numpy array of data with dimensions (lead_time, lat, lon, variable, member*time)
   _data_agg      A 3D numpy array of data with dimensions (lead_time, variable, member*time)
   """
   def __init__(self):
      self._debug = False
      """ External state (such as the month of the year) """
      self._ext_state = None
      self.aggregator = wxgen.aggregator.Mean()
      self._data_agg_cache = None

   def info(self):
      print "Database information:"
      print "  Length of segments: %d" % self.length
      print "  Number of segments: %d" % self.num
      print "  Number of variables: %d" % len(self.variables)

   def get(self, i):
      """ Get the i'th trajectory in the database """
      indices = np.zeros([self.length, 2], int)
      indices[:,0] = i
      indices[:,1] = np.arange(0, self.length)
      return wxgen.trajectory.Trajectory(indices, self)

   def get_truth(self):
      """ Concatenate the initialization state of all trajectories """
      times = np.arange(np.min(self.inittimes), np.max(self.inittimes), 86400)
      indices = np.nan*np.zeros([len(times), 2], int)
      for i in range(0, len(times)):
         time = times[i]
         I = np.where(time >= self.inittimes)[0]
         if len(I) == 0:
            wxgen.error("Internal error")
         inittime = self.inittimes[I[-1]]
         lt = (time - inittime)/86400
         if lt < self.length:
            indices[i,0] = I[-1]
            indices[i,1] = lt
      return wxgen.trajectory.Trajectory(indices, self)

   def get_random(self, target_state, metric, ext_state=None):
      """
      Returns a random segment from the database that is weighted
      by the scores computed by metric.

      Arguments:
      target_state   A numpy array (length V)
      metric         Of type wxgen.metric.Metric
      ext_state      External state

      Returns:
      trajectory     Of type wxgen.trajectory.Trajectory
      """
      weights = metric.compute(target_state, self._data_agg[0,:,:])
      Ivalid = np.where(np.isnan(weights) == 0)[0]
      if ext_state is not None and self._ext_state is not None:
         Iext_state = np.where(self._ext_state[Ivalid] == ext_state)[0]
         if len(Iext_state) == 0:
            wxgen.util.error("Cannot find a segment with  external state = %s" % str(ext_state))
         Ivalid = Ivalid[Iext_state]

      weights_v = weights[Ivalid]

      # Flip the metric if it is negative oriented
      if metric._orientation == -1:
         I0 = np.where(weights_v < 1e-3)[0]
         I1 = np.where(weights_v >= 1e-3)[0]
         # Ensure we do not get too high weights
         weights_v[I1] = 1.0/weights_v[I1]
         weights_v[I0] = 1e3

      I_v = wxgen.util.random_weighted(weights_v)
      I = Ivalid[I_v]

      # Do a weighted random choice of the weights
      if self._debug:
         print "Data: ", self._data_agg[0,:,I]
         print "Weight: ", weights_v[I_v]
         print "Max weight: ", np.max(weights_v)
      return self.get(I)

   def get_sequence(self, indices):
      """ Returns a gridded sequence of states
      
      Arguments:
      indices     A numpy array of integers

      Returns:
      data        A 4D array (T, X, Y, V)
      """
      pass

   @property
   def _data_agg(self):
      if self._data_agg_cache is None:
         self._data_agg_cache = np.zeros([self.length, len(self.variables), self.num], float)
         self._data_agg_cache = self.aggregator(self._data)
      return self._data_agg_cache

   @property
   def X(self):
      return self._data.shape[1]

   @property
   def Y(self):
      return self._data.shape[2]

   def index2(self):
      pass


class Random(Database):
   """
   Trajectories are random Gaussian walks, with constant variance over time. There is no covariance
   between different forecast variables.
   """
   def __init__(self, N, T, V, variance=1):
      Database.__init__(self)
      self.num = N
      self.length = T
      if V == None:
         V = 1
      self._V = V
      self._variance = variance
      self._data = np.zeros([T, 1, 1, V, N], float)

      # Ensure that the signal has a constant variance over time
      scale = 1./np.sqrt(np.linspace(1, T, T))

      for v in range(0, self._V):
         self._data[:,0,0,v,:]  = np.transpose(np.resize(scale, [N, T])) * np.cumsum(np.random.randn(T, N)*np.sqrt(self._variance), axis=0)

      self.variables = [wxgen.variable.Variable(str(i)) for i in range(0, self._V)]
      self.lats = [0]
      self.lons = [0]


class Netcdf(Database):
   """
   Segments stored in a netcdf database

   Should have the following format:
      dims: date, leadtime, member
      vars: date(date), leadtime(leadtime)
            variable_name(date, leadtime, member)
   where variable_name is one or more names of weather variables

   Internal
   _data0         A 6D numpy array of data with dimensions (time, lead_time, ensemble_member, lat, lon, variable)
   """
   _debug0 = True
   def __init__(self, filename, V=None):
      """
      filename    Load data from this file
      V           Only use the first V variables in the database
      """
      Database.__init__(self)
      self._file = netCDF4.Dataset(filename)

      self._initname = "forecast_reference_time"

      # Set dimensions
      var_names = [name for name in self._file.variables if name not in ["lat", "lon", "ensemble_member", "time", "dummy", "longitude_latitude", "forecast_reference_time"]]
      self.variables = list()
      for var_name in var_names:
         units = None
         label = None
         if hasattr(self._file.variables[var_name], "units"):
            units = self._file.variables[var_name].units
         if hasattr(self._file.variables[var_name], "standard_name"):
            label = self._file.variables[var_name].standard_name
         var = wxgen.variable.Variable(var_name, units, label)
         self.variables.append(var)
      if V is not None:
         self.variables = self.variables[0:V]

      # Load data
      self.length = self._file.dimensions["time"].size
      self._members = self._file.dimensions["ensemble_member"].size
      V = len(self.variables)
      M = self._members
      D = self._file.dimensions["forecast_reference_time"].size
      T = self.length
      if "lon" in self._file.dimensions:
         is_spatial = True
         X = self._file.dimensions["lon"].size
         Y = self._file.dimensions["lat"].size
         self.lats = self._copy(self._file.variables["lat"])
         self.lons = self._copy(self._file.variables["lon"])
      else:
         is_spatial = False
         X = 1
         Y = 1
         self.lats = [0]
         self.lons = [0]
      self.num = M * D
      print "Allocating %.2f GB" % (T*Y*X*V*M*D*4.0/1024/1024/1024)
      self._data = np.nan*np.zeros([T, Y, X, V, M*D], float)
      self._ext_state = np.zeros(self.num, float)
      assert(self._ext_state.shape[0] == self.num)

      for v in range(0, V):
         var = self.variables[v]
         print var.name
         temp = self._copy(self._file.variables[var.name]) # dims: D, T, M, X, Y)

         # Quality control
         if var.name == "precipitation_amount":
            temp[temp < 0] = np.nan
         index = 0
         for d in range(0, D):
            for m in range(0, M):
               if is_spatial:
                  self._data[:,:,:,v,index] = temp[d, :, m, :, :]
               else:
                  self._data[:,:,:,v,index] = np.reshape(temp[d, :, m], [T,Y,X])
               index = index + 1

      # If one or more values are missing for a member, set all values to nan
      for e in range(0, M*D):
         if np.sum(np.isnan(self._data[:,:,:,:,e])) > 0:
            self._data[:,:,:,:,e] = np.nan
            # print "Removing %d" % e

      times = self._file.variables[self._initname]
      day_of_year = np.zeros(times.shape)
      for d in range(0, times.shape[0]):
         day_of_year[d] = int(datetime.datetime.fromtimestamp(times[d]).strftime('%j'))
      month_of_year = day_of_year.astype(int)/ 30
      self.inittimes = np.repeat(times, self._members)
      self._ext_state = np.repeat(month_of_year, self._members)
      if self._debug0:
         print self._ext_state

      self._file.close()

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
      self.num = N
      self.length = T
      self._R = R
      self._S = S
      self._B = B
      self._dt = dt
      self._V = 3  # Number of variables
      self._data = np.zeros([T, 1, 1, self._V, N], float)
      self.lats = [0]
      self.lons = [0]
      self._initial_state = [-10, -10, 25]
      self._std_initial_state = 0.1  # Standard deviation of initial condition error

      # Initialize
      for v in range(0, self._V):
         self._data[0,0,0,v,:] = self._initial_state[v] + np.random.randn(N) * self._std_initial_state

      TT = int(1 / self._dt)/10
      # Iterate
      for t in range(1, T):
         x0 = copy.deepcopy(self._data[t-1,0,0,0,:])
         y0 = copy.deepcopy(self._data[t-1,0,0,1,:])
         z0 = copy.deepcopy(self._data[t-1,0,0,2,:])
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

         self._data[t,0,0,0,:] = x0
         self._data[t,0,0,1,:] = y0
         self._data[t,0,0,2,:] = z0

      self.variables = [wxgen.variable.Variable(i) for i in ["X", "Y", "Z"]]
