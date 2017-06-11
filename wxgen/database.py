import numpy as np
import wxgen.util
import pywt
import os
import datetime
import wxgen.variable
import wxgen.climate_model
import copy
import netCDF4
import time as timing


class Database(object):
   """
   Abstract class which stores weather segments (short time-series)

   The database is stored in self._data_cache, which is a dictionary with variable key and arrays of
   dimensions (T, X, Y, N) as values, where T is the segment length, X and Y are geographical
   dimensions, and N is the number of segments.

   Attributes:
      variables: A list of wxgen.variable.Variable
      length (int): Length of each trajectory (in number of days)
      num (int): Number of trajectories
      lats (np.array): 2D grid of latitudes
      lons (np.array): 2D grid of longitudes
      X (int): Number of X-axis points
      Y (int): Number of Y-axis points
      inittimes (np.array): array with initialization time corresponding to each member
      climate_states (np.array): array of climate states (one for each ensemble member)
      name: Name of this database (e.g. filename)

   Subclass requirements (these must be filled out):
      self.name
      self.length
      self.num
      self.variables
      self.inittimes
      self.lats
      self.lons
      _load(self, variable)

   """
   def __init__(self, model=None):
      self._data_matching_cache = None
      self.wavelet_levels = 0
      self.mem = None
      if model is None:
         self.model = wxgen.climate_model.Bin(10)
      else:
         self.model = model

      # Cache variables to store data returned by @property functions
      self._data_cache = dict()
      self._data_agg_cache = None
      self._climate_states_cache = None

   def info(self):
      print "Database information:"
      print "  Length of segments: %d" % self.length
      print "  Number of segments: %d" % self.num
      print "  Number of variables: %d" % len(self.variables)

   def load(self, variable):
      """
      Loads the variable
      """
      t = timing.time()
      if variable not in self._data_cache:
         wxgen.util.debug("Cache miss variable '%s'" % variable.name)

         if len(self._data_cache) > 0:
            # Check if we need to remove data from cache
            akey = self._data_cache.keys()[0]
            bytes_per_value = 4
            size_per_key = np.product(self._data_cache[akey].shape) * bytes_per_value
            next_size = float(len(self._data_cache) + 1) * size_per_key
            next_size_gb = next_size / 1e9
            if self.mem is not None and next_size_gb > self.mem:
               # remove from cache
               I = np.random.randint(len(self._data_cache))
               rmkey = self._data_cache.keys()[I]
               self._data_cache.pop(rmkey)
               wxgen.util.warning("Cache full (%2.1fGB): Removing member '%s' from cache" % (next_size_gb, rmkey.name))

         # Get data from subclass
         data = self._load(variable)
         self._data_cache[variable] = data
         e = timing.time()
         # print "Timing: %g" % (e - t)
      else:
         wxgen.util.debug("Cache hit '%s'" % variable.name)

      return self._data_cache[variable]

   def _load(self, variable):
      raise NotImplementedError()

   def get(self, i):
      """ Get the i'th trajectory in the database """
      indices = np.zeros([self.length, 2], int)
      indices[:, 0] = i
      indices[:, 1] = np.arange(0, self.length)
      assert(np.sum(np.isnan(indices)) == 0)
      return wxgen.trajectory.Trajectory(indices)

   def get_truth(self, start_date=None, end_date=None):
      """ Concatenate the initialization state of all trajectories """
      start = np.min(self.inittimes) if start_date is None else wxgen.util.date_to_unixtime(start_date)
      end = np.max(self.inittimes) if end_date is None else wxgen.util.date_to_unixtime(end_date)
      times = np.arange(start, end, 86400)
      indices = -1*np.ones([len(times), 2], int)
      wxgen.util.debug("Start: %d End: %d" % (start, end))
      for i in range(0, len(times)):
         time = times[i]
         I = np.where(time >= self.inittimes)[0]
         if len(I) == 0:
            wxgen.util.error("There are no inittimes available before %d. The earliest is %d." %
                  (wxgen.util.unixtime_to_date(time), wxgen.util.unixtime_to_date(
                     np.min(self.inittimes))))
         inittime = np.max(self.inittimes[I])
         lt = int((time - inittime)/86400)
         if lt < self.length:
            indices[i, 0] = np.where(self.inittimes == inittime)[0][0]
            indices[i, 1] = lt
         else:
            wxgen.util.warning("Did not find an index for %d = %d. Using the previous state." % (time, wxgen.util.unixtime_to_date(time)))
            assert(i > 0)
            indices[i, 0] = indices[i-1, 0]
            indices[i, 1] = indices[i-1, 1]

      return wxgen.trajectory.Trajectory(indices)

   def extract(self, trajectory):
      """
      Extract a trajectory of large-scale aggregated values from the database

      Arguments:
         trajectory (wxgen.trajectory.Trajectory): Trajectory to extract

      Returns:
         np.array: A 2D array (Time, variable) sequence of values
      """
      T = trajectory.indices.shape[0]
      V = len(self.variables)
      values = np.nan*np.zeros([T, V], float)
      for i in range(0, trajectory.indices.shape[0]):
         if trajectory.indices[i, 1] >= 0:
            values[i, :] = self._data_agg[trajectory.indices[i, 1], :, trajectory.indices[i, 0]]
      return values

   def extract_grid(self, trajectory, variable):
      """
      Extract a trajectory of large-scale values from the database

      Arguments:
         trajectory (wxgen.trajectory.Trajectory): Trajectory to extract
         variable (wxgen.variable.Variable): Variable to extract

      Returns:
         np.array: A 4D array (Time, X, Y) sequence of values
      """
      T = trajectory.indices.shape[0]
      V = len(self.variables)
      X = self.X
      Y = self.Y
      temp = self.load(variable)
      values = np.nan*np.zeros([T, Y, X], float)
      # Loop over member, lead-time indices
      st = timing.time()
      for i in range(0, trajectory.indices.shape[0]):
         m = trajectory.indices[i, 0]
         t = trajectory.indices[i, 1]
         # print i, m, t
         assert(not np.isnan(m))
         assert(not np.isnan(t))
         if t >= 0:
            values[i, :, :] = temp[t, :, :, m]
      e = timing.time()
      # print "Q Timing: %g" % (e - st)
      return values

   def extract_matching(self, trajectory):
      """
      Extract a trajectory of values used to match states from the database

      Arguments:
         trajectory (wxgen.trajectory.Trajectory): Trajectory to extract

      Returns:
         np.array: A 2D array (Time, variable) sequence of values
      """
      T = trajectory.indices.shape[0]
      values = np.nan*np.zeros([T, self._data_matching.shape[1]], float)
      for i in range(0, trajectory.indices.shape[0]):
         if trajectory.indices[i, 1] >= 0:
            values[i, :] = self._data_matching[trajectory.indices[i, 1], :, trajectory.indices[i, 0]]
      return values

   def get_wavelet_size(self):
      if self.wavelet_levels == 0:
         return 1, 1
      NX = int(np.ceil(float(self.Y)/2**self.wavelet_levels))
      NY = int(np.ceil(float(self.X)/2**self.wavelet_levels))
      return NX, NY

   @property
   def _data_agg(self):
      # _data_agg (np.array): A 3D array of data with dimensions (lead_time, variable, member*time)
      if self._data_agg_cache is None:
         self._data_agg_cache = np.zeros([self.length, self.V, self.num])
         for v in range(len(self.variables)):
            variable = self.variables[v]
            data = self.load(variable)
            self._data_agg_cache[:, v, :] = np.mean(np.mean(data, axis=2), axis=1)

      return self._data_agg_cache

   @property
   def climate_states(self):
      if self._climate_states_cache is None:
         self._climate_states_cache = self.model.get(self.inittimes)
      return self._climate_states_cache

   @property
   def V(self):
      return len(self.variables)

   @property
   def X(self):
      return self.lats.shape[1]

   @property
   def Y(self):
      return self.lats.shape[0]

   @property
   def _data_matching(self):
      """
      Get the data used to match states. This may be different than data_agg, since it can include
      wavelet information.

      Returns:
         np.array: Dimensions self.length, variable, self.num
      """
      if self._data_matching_cache is None:
         if self.wavelet_levels == 0:
            self._data_matching_cache = self._data_agg
         else:
            """
            Decompose all gridded fields into wavelet components. Do this separately for each
            variable, leadtime, and ensemble member and store the values in
            self._data_matching_cache.
            """
            wxgen.util.debug("Computing wavelets")
            s = timing.time()
            NX, NY = self.get_wavelet_size()
            N = int(NX * NY)
            self._data_matching_cache = np.zeros([self.length, self.V*N, self.num])
            for v in range(self.V):
               data = self.load(self.variables[v])
               # Compute wavelet decomposition on axes 1 2
               dec = pywt.wavedec2(data, 'haar', level=self.wavelet_levels, axes=(1,2))[0]
               dec = dec / 2**self.wavelet_levels

               # Collapse axes 1 and 2
               dec = np.reshape(dec, [self.length, N, self.num])

               # Write data into the right location in the cache
               I = range(v*N, (v+1)*N)
               self._data_matching_cache[:, I, :] = dec

            e = timing.time()
            wxgen.util.debug("Wavelet time: %f" % (e - s))

      return self._data_matching_cache


class Netcdf(Database):
   """
   Segments stored in a netcdf database. This can either be an aggregated or a gridded database.

   Aggregated database has the following format:
      dims: forecast_reference_time, time, ensemble_member
      vars: forecast_reference_time(forecast_reference_time), time(time)
            variable_name(forecast_reference_time, time, ensemble_member)

   Gridded database has the following format:
      dims: forecast_reference_time, lat(itude), lon(gitude), time, ensemble_member
      vars: forecast_reference_time(forecast_reference_time), time(time)
            variable_name(forecast_reference_time, time, member, latitude, longitude)

   where variable_name is one or more names of weather variables. Forecast_reference_time is
   optional, i.e. both the variable and dimension could be missing.
   """
   def __init__(self, filename, vars=None, model=None, mem=None):
      """
      Arguments:
         filename (str): Load data from this file
         vars (list): List of indices for which variables to use
      """
      Database.__init__(self, model)
      if not os.path.isfile(filename):
         wxgen.util.error("File '%s' does not exist" % filename)

      self.name = filename[filename.rfind('/') + 1:]

      self.mem = mem
      self._file = netCDF4.Dataset(filename)

      # Set dimensions
      var_names = [name for name in self._file.variables if name not in ["lat", "lon", "latitude", "longitude", "x", "y", "ensemble_member", "time", "dummy", "longitude_latitude", "forecast_reference_time", "projection_regular_ll"]]
      if vars is None:
         vars = range(len(var_names))

      if vars is not None and max(vars) >= len(var_names):
         wxgen.util.error("Index in --vars (%d) is >= number of variables (%d)" % (max(vars), len(var_names)))

      self.variables = list()
      for i in vars:
         var_name = var_names[i]
         units = None
         label = None
         dims = self._file.variables[var_name].dimensions
         if "time" in dims:
            if hasattr(self._file.variables[var_name], "units"):
               units = self._file.variables[var_name].units
            if hasattr(self._file.variables[var_name], "standard_name"):
               label = self._file.variables[var_name].standard_name
            var = wxgen.variable.Variable(var_name, units, label)
            self.variables.append(var)
            wxgen.util.debug("Using variable '%s'" % var_name)

      # Load data
      self.length = self._file.dimensions["time"].size

      self.has_frt = True
      if "forecast_reference_time" in self._file.dimensions:
         times = self._file.variables["forecast_reference_time"][:]
         self.ens = self._file.dimensions["ensemble_member"].size
         self.num = self._file.dimensions["forecast_reference_time"].size * self.ens
         self.inittimes = np.repeat(times, self.ens)
      else:
         self.has_frt = False
         times = np.array([self._file.variables["forecast_reference_time"][:]])
         self.num = self._file.dimensions["ensemble_member"].size
         self.inittimes = times
         self.ens = self.num

      self._Itimes = np.where(np.isnan(times) == 0)[0]
      times = times[self._Itimes]

      # Read lat/lon dimensions
      self.is_spatial = True
      if "lon" in self._file.dimensions:
         X = self._file.dimensions["lon"].size
         Y = self._file.dimensions["lat"].size
      elif "longitude" in self._file.dimensions:
         X = self._file.dimensions["longitude"].size
         Y = self._file.dimensions["latitude"].size
      elif "x" in self._file.dimensions:
         X = self._file.dimensions["x"].size
         Y = self._file.dimensions["y"].size
      else:
         self.is_spatial = False
         X = 1
         Y = 1
         self.lats = np.zeros([1, 1])
         self.lons = np.zeros([1, 1])

      # Read lat/lon variables
      if self.is_spatial:
         if "lat" in self._file.variables:
            self.lats = self._copy(self._file.variables["lat"])
            self.lons = self._copy(self._file.variables["lon"])
         elif "latitude" in self._file.variables:
            self.lats = self._copy(self._file.variables["latitude"])
            self.lons = self._copy(self._file.variables["longitude"])
         if len(self.lats.shape) == 1 and len(self.lons.shape) == 1:
            wxgen.util.debug("Meshing latitudes and longitudes")
            self.lons, self.lats = np.meshgrid(self.lons, self.lats)

   def _load(self, variable):
      temp = self._file.variables[variable.name][:]

      data = np.nan*np.zeros([self.length, self.Y, self.X, self.num], np.float32)

      index = 0
      for d in range(len(self._Itimes)):
         for m in range(0, self.ens):
            if self.is_spatial:
               if self.has_frt:
                  data[:, :, :, index] = temp[self._Itimes[d], :, m, :, :]
               else:
                  data[:, :, :, index] = temp[:, m, :, :]
            else:
               if self.has_frt:
                  data[:, :, :, index] = np.reshape(temp[self._Itimes[d], :, m], [self.length, self.Y, self.X])
               else:
                  data[:, :, :, index] = np.reshape(temp[:, m], [self.length, self.Y, self.X])
            # If one or more values are missing for a member, set all values to nan
            NM = np.sum(np.isnan(data[:, :, :, index]))
            if NM > 0:
               data[:, :, index] = np.nan
               wxgen.util.debug("Removing member %d because of %d missing values" % (index, NM))
            index = index + 1
      # Quality control
      if variable.name == "precipitation_amount":
         data[data < 0] = 0

      return data

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
   def __init__(self, num=1000, length=100, R=28, S=10, B=2.6667, dt=0.0001, model=wxgen.climate_model.Zero()):
      Database.__init__(self, model)
      self.name = "Lorenz(%d,%d,%f,%f,%f)" % (num, length, R, S, B)
      self.length = length
      self.num = num
      self._R = R
      self._S = S
      self._B = B
      self._dt = dt
      self.lats = np.zeros([1, 1])
      self.lons = np.zeros([1, 1])
      var_x = wxgen.variable.Variable("X")
      var_y = wxgen.variable.Variable("Y")
      var_z = wxgen.variable.Variable("Z")
      self.variables = [var_x, var_y, var_z]
      self._initial_state = {var_x: -10, var_y: -10, var_z: -25}
      self._std_initial_state = 10  # Standard deviation of initial condition error
      self.inittimes = np.ones(self.num) + wxgen.util.date_to_unixtime(20150101)

      # Initialize
      self._data = dict()
      for var in self.variables:
         self._data[var] = np.zeros([self.length, self.X, self.Y, self.num])
         self._data[var][0, 0, 0, :] = self._initial_state[var] + np.random.randn(self.num) * self._std_initial_state

      TT = int(1 / self._dt)/10
      # Iterate
      for t in range(1, self.length):
         x0 = copy.deepcopy(self._data[var_x][t-1, 0, 0, :])
         y0 = copy.deepcopy(self._data[var_y][t-1, 0, 0, :])
         z0 = copy.deepcopy(self._data[var_z][t-1, 0, 0, :])
         # Iterate from t-1 to t
         for tt in range(0, TT):
            x1 = x0 + self._dt * self._S * (y0 - x0)
            y1 = y0 + self._dt * (self._R * x0 - y0 - (x0 * z0))
            z1 = z0 + self._dt * (x0*y0 - self._B*z0)

            x2 = x1 + self._dt * self._S * (y1 - x1)
            y2 = y1 + self._dt * (self._R * x1 - y1 - x1*z1)
            z2 = z1 + self._dt * (x1*y1 - self._B*z1)

            x0 = 0.5 * (x2 + x0)
            y0 = 0.5 * (y2 + y0)
            z0 = 0.5 * (z2 + z0)

         self._data[var_x][t, 0, 0, :] = x0
         self._data[var_y][t, 0, 0, :] = y0
         self._data[var_z][t, 0, 0, :] = z0

   def _load(self, variable):
      return self._data[variable]


class Random(Database):
   """
   A minimal database to check that the subclass requirements specified are correct

   Trajectories are random Gaussian walks, with constant variance over time. There is no covariance
   between different forecast variables.
   """
   def __init__(self, num=365, length=10, num_vars=1, model=None):
      Database.__init__(self, model)
      self.name = "Test"
      self.length = length
      self.num = num
      self.variables = [wxgen.variable.Variable("temperature")] * num_vars
      self.inittimes = np.array(range(self.num)) * 86400 + wxgen.util.date_to_unixtime(20150101)
      X = 10
      Y = 8
      self.lats, self.lons = np.meshgrid(range(50, 50 + Y), range(X))

      self._variance = 1

   def _load(self, variable):
      values = np.random.randn(self.length, self.Y, self.X, self.num) * np.sqrt(self._variance)
      values = np.cumsum(values, axis=0)
      for i in range(self.length):
         # Ensure that the signal has a constant variance over time
         scale = 1./np.sqrt(1 + i)
         values[i, :, :, :] = values[i, :, :, :] * scale
      return values
