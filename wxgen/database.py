import numpy as np
import wxgen.util
import pywt
import os
import sys
import datetime
import wxgen.variable
import wxgen.climate_model
import copy
import netCDF4
import time as timing
import wxgen.config


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
      altitudes (np.array): 2D grid of altitudes
      X (int): Number of X-axis points
      Y (int): Number of Y-axis points
      inittimes (np.array): array with initialization time corresponding to each member
      leadtimes (np.array): array with leadtimes in seconds
      climate_states (np.array): array of climate states (one for each ensemble member)
      name: Name of this database (e.g. filename)

   Subclass requirements (these must be filled out):
      self.name
      self.length
      self.num
      self.variables
      self.inittimes
      self.leadtimes
      self.lats
      self.lons
      self.altitudes
      self.timestep (Number of seconds between timesteps, e.g. 86400 or 6 * 3600)
      _load(self, variable)

   """
   def __init__(self, model=None):
      self._data_matching_cache = None
      self.spatial_decomposition = 0
      self.join_config = None
      self.mem = None
      if model is None:
         self.model = wxgen.climate_model.Bin(10)
      else:
         self.model = model

      # Cache variables to store data returned by @property functions
      self._data_cache = dict()
      self._data_agg_cache = None
      self._climate_states_cache = None
      self._label = None

   @property
   def label(self):
      """
      Use this to label diagrams. The label can overwritten by the application if needed, by default
      the name of the database is used.
      """
      if self._label is not None:
         return self._label
      else:
         return self.name

   @label.setter
   def label(self, value):
      self._label = value

   def info(self):
      print("Database information:")
      print("  Length of segments: %d" % self.length)
      print("  Number of segments: %d" % self.num)
      print("  Number of variables: %d" % len(self.variables))
      unique_states = np.sort(np.unique(self.climate_states))
      print("  Climate states: " + ', '.join([str(s) for s in unique_states]))

   def load(self, variable):
      """
      Loads the variable
      """
      t = timing.time()
      if variable not in self._data_cache:
         wxgen.util.debug("Cache miss variable '%s'" % variable.name)

         if len(self._data_cache) > 0:
            # Check if we need to remove data from cache
            akey = list(self._data_cache.keys())[0]
            bytes_per_value = 4
            size_per_key = np.product(self._data_cache[akey].shape) * bytes_per_value
            next_size = float(len(self._data_cache) + 1) * size_per_key
            next_size_gb = next_size / 1e9
            if self.mem is not None and next_size_gb > self.mem:
               # remove from cache
               I = np.random.randint(len(self._data_cache))
               rmkey = list(self._data_cache.keys())[I]
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

   def get_truth(self, start_date=None, end_date=None, start_time_of_day=None, which_leadtime=0):
      """
      Concatenate the initialization state of all trajectories

      Arguments:
         start_date (int): Starting date of scenario (YYYYMMDD)
         end_date (int): Starting date of scenario (YYYYMMDD)
         which_leadtime (int): Which lead time should be used to create the truth?

      Returns:
         trajectory.Trajectory: Truth trajectory
      """
      start = np.min(self.inittimes) if start_date is None else wxgen.util.date_to_unixtime(start_date)
      end = np.max(self.inittimes) if end_date is None else wxgen.util.date_to_unixtime(end_date)
      times = np.arange(start, end, self.timestep)
      indices = -1*np.ones([len(times), 2], int)
      wxgen.util.debug("Start: %d End: %d" % (start, end))
      for i in range(0, len(times)):
         time = times[i]
         I = np.where(time >= self.inittimes)[0]
         if len(I) == 0:
            wxgen.util.error("There are no inittimes available before %d. The earliest is %d." %
                  (wxgen.util.unixtime_to_date(time), wxgen.util.unixtime_to_date(
                     np.min(self.inittimes))))
         """ Find the inittime and leadtime that is closest to the desired date and 'which_leadtime' """
         curr_times = self.inittimes[I] + self.timestep * which_leadtime
         Ibest = np.argmin(np.abs(curr_times - time))
         inittime = self.inittimes[I][Ibest]
         lt = int((time - inittime)/self.timestep)
         if lt < self.length:
            indices[i, 0] = np.where(self.inittimes == inittime)[0][0]
            indices[i, 1] = lt
         else:
            wxgen.util.warning("Did not find an index for %d = %d. Using the previous state." % (time, wxgen.util.unixtime_to_date(time)))
            assert(i > 0)
            indices[i, 0] = indices[i-1, 0]
            indices[i, 1] = indices[i-1, 1]

      if start_time_of_day is not None:
         # Crop away the first few timesteps to make the trajectory start at the right time of day
         Iindices = np.where(self.leadtimes[indices[:, 1]] % 86400 == start_time_of_day)[0][0]
         indices = indices[Iindices:, :]
         # TODO: Deal with the fact that the cropping makes it so that the trajectory is too short

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
      values = np.nan*np.zeros([T, self._data_matching.shape[1]], 'float32')
      for i in range(0, trajectory.indices.shape[0]):
         if trajectory.indices[i, 1] >= 0:
            values[i, :] = self._data_matching[trajectory.indices[i, 1], :, trajectory.indices[i, 0]]
      return values

   def get_wavelet_size(self):
      if self.spatial_decomposition == 0:
         return 1, 1
      elif self.spatial_decomposition == 'all':
         return self.X, self.Y
      else:
         NX = int(np.ceil(float(self.Y)/2**self.spatial_decomposition))
         NY = int(np.ceil(float(self.X)/2**self.spatial_decomposition))
         return NX, NY

   @property
   def _data_agg(self):
      # _data_agg (np.array): A 3D array of data with dimensions (lead_time, variable, member*time)
      if self._data_agg_cache is None:
         self._data_agg_cache = np.zeros([self.length, self.V, self.num], 'float32')
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
         if self.join_config is not None:
            config = wxgen.config.Config(self.join_config)

            # Preload all required variables
            Ivars = set([int(point['variable']) for point in config.points])
            data = dict()
            for Ivar in Ivars:
               data[Ivar] = self.load(self.variables[Ivar])

            self._data_matching_cache = np.zeros([self.length, len(config.points), self.num], 'float32')
            for p, point in enumerate(config.points):
               # Find nearest neighbour
               dist = (self.lats - point['lat']) ** 2 + (self.lons - point['lon']) ** 2
               indices = np.unravel_index(dist.argmin(), dist.shape)
               Ivar = point['variable']
               weight = point['weight']
               self._data_matching_cache[:, p, :] = data[Ivar][:, indices[0], indices[1], :] * weight
         elif self.spatial_decomposition == 0:
            self._data_matching_cache = self._data_agg
         elif self.spatial_decomposition == 'all':
            N = int(self.X * self.Y)
            self._data_matching_cache = np.zeros([self.length, self.V*N, self.num], 'float32')
            for v in range(self.V):
               data = self.load(self.variables[v])
               data = data[:, :, 0, :]
               # Write data into the right location in the cache
               I = range(v*N, (v+1)*N)
               self._data_matching_cache[:, I, :] = data
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
               dec = pywt.wavedec2(data, 'haar', level=self.spatial_decomposition, axes=(1, 2))[0]
               dec = dec / 2**self.spatial_decomposition

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

   Concatenated database has the following format:
      dims: time, grid_point, lead_time, ensemble_member
      vars: time(time), lead_time(lead_time)
            variable_name(time, lead_time, ensemble_member, grid_point)

   Gridded database has the following format:
      dims: time, lat(itude), lon(gitude), lead_time, ensemble_member
      vars: time(time), lead_time(lead_time)
            variable_name(time, lead_time, ensemble_member, latitude, longitude)

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
      var_names = [name for name in self._file.variables if name not in ["lat", "lon", "latitude",
      "longitude", "altitude", "x", "y", "ensemble_member", "lead_time", "dummy",
      "longitude_latitude", "time", "projection_regular_ll", "segment_lead_time", "segment_member", "segment_time"]]
      if vars is None:
         vars = range(len(var_names))

      if vars is not None and max(vars) >= len(var_names):
         wxgen.util.error("Input has only %d variables. Variable %d is outside range." % (len(var_names), max(vars)))

      if "lead_time" in self._file.dimensions and "time" in self._file.dimensions:
         lead_time_dim = "lead_time"
         time_dim = "time"
         self.single_run = False
      elif "time" in self._file.dimensions:
         lead_time_dim = "time"
         time_dim = None
         self.single_run = True
      else:
         wxgen.util.error("Cannot read file. Missing 'time' and/or 'lead_time' dimensions")
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
      self.length = len(self._file.dimensions[lead_time_dim])
      timevar = self._file.variables[lead_time_dim]

      # Assume dataset has a daily timestep, except if it is possible to deduce otherwise
      self.timestep = 86400
      leadtimes = wxgen.util.clean(self._file.variables[lead_time_dim][:])
      if hasattr(timevar, "units"):
         if len(timevar.units) >= 7 and timevar.units[0:7] == "seconds":
            pass
         elif timevar.units == "days":
            leadtimes *= 86400
         elif timevar.units == "hours":
            leadtimes *= 3600
         else:
            wxgen.error("Cannot parse time units")
      self.timestep = leadtimes[1] - leadtimes[0]
      self.leadtimes = leadtimes
      if np.isnan(self.timestep):
         wxgen.util.error("Cannot determine timestep from database")

      self.has_single_spatial_dim = False
      if self.single_run:
         if "forecast_reference_time" in self._file.variables:
            times = np.array([self._file.variables["forecast_reference_time"][:]])
         else:
            times = np.array([self._file.variables["time"][0]])
         self.num = len(self._file.dimensions["ensemble_member"])
         self.inittimes = times
         self.ens = self.num
      else:
         times = self._file.variables["time"][:]
         self.ens = len(self._file.dimensions["ensemble_member"])
         self.num = len(self._file.dimensions["time"]) * self.ens
         self.inittimes = np.repeat(times, self.ens)

      self._Itimes = np.where(np.isnan(times) == 0)[0]
      times = times[self._Itimes]

      # Read lat/lon dimensions
      if "lon" in self._file.dimensions:
         X = len(self._file.dimensions["lon"])
         Y = len(self._file.dimensions["lat"])
      elif "longitude" in self._file.dimensions:
         X = len(self._file.dimensions["longitude"])
         Y = len(self._file.dimensions["latitude"])
      elif "x" in self._file.dimensions:
         X = len(self._file.dimensions["x"])
         Y = len(self._file.dimensions["y"])
      elif "grid_point" in self._file.dimensions:
         # Store concatinated spatial dimension in X dimension
         X = len(self._file.dimensions["grid_point"])
         Y = 1
         self.has_single_spatial_dim = True
      else:
         raise NotImplementedError

      """
      Read lat/lon/elev variables
      """
      if self.has_single_spatial_dim:
         self.lats = self._copy(self._file.variables["latitude"])
         self.lons = self._copy(self._file.variables["longitude"])
         self.lats = np.expand_dims(self.lats, 1)
         self.lons = np.expand_dims(self.lons, 1)
      else:
         if "lat" in self._file.variables:
            self.lats = self._copy(self._file.variables["lat"])
            self.lons = self._copy(self._file.variables["lon"])
         elif "latitude" in self._file.variables:
            self.lats = self._copy(self._file.variables["latitude"])
            self.lons = self._copy(self._file.variables["longitude"])
         if len(self.lats.shape) == 1 and len(self.lons.shape) == 1:
            wxgen.util.debug("Meshing latitudes and longitudes")
            self.lons, self.lats = np.meshgrid(self.lons, self.lats)

      """
      Read altitude information
      """
      if "altitude" in self._file.variables:
         self.altitudes = self._copy(self._file.variables["altitude"])
         if self.altitudes.shape != self.lats.shape:
            # Try to remove singleton dimensions
            self.altitudes = np.squeeze(self.altitudes)
      elif "surface_geopotential" in self._file.variables:
         self.altitudes = np.squeeze(self._copy(self._file.variables["surface_geopotential"]/9.81))
         if self.altitudes.shape != self.lats.shape:
            self.altitudes = np.squeeze(self.altitudes)
      else:
         self.altitudes = np.nan * self.lats
      if self.has_single_spatial_dim:
         self.altitudes = np.expand_dims(self.altitudes, 1)
      if self.altitudes.shape != self.lons.shape:
         wxgen.util.error("Altitude dimensions do not match those of lat/lon")

   def _load(self, variable):
      if variable.name not in self._file.variables:
         wxgen.util.error("Variable '%s' does not exist in file '%s'" % (variable.name, self.name))
      wxgen.util.debug("Allocating %g GB for '%s'" % (np.product(self._file.variables[variable.name].shape) * 4.0 / 1e9, variable.name))
      temp = self._file.variables[variable.name][:]

      data = np.nan * np.zeros([self.length, self.Y, self.X, self.num], 'float32')

      index = 0
      for d in range(len(self._Itimes)):
         for m in range(0, self.ens):
            if self.has_single_spatial_dim:
               if self.single_run:
                  data[:, :, 0, index] = temp[:, m, :]
               else:
                  data[:, :, 0, index] = temp[self._Itimes[d], :, m, :]
            else:
               if self.single_run:
                  data[:, :, :, index] = temp[:, m, :, :]
               else:
                  data[:, :, :, index] = temp[self._Itimes[d], :, m, :, :]
            # If one or more values are missing for a member, set all values to nan
            NM = np.sum(np.isnan(data[:, :, :, index]))
            if NM > 0:
               data[:, :, index] = np.nan
               wxgen.util.debug("Removing member %d because of %d missing values" % (index, NM))
            index = index + 1
      # Quality control
      if variable.name == "precipitation_amount":
         data[data < 0] = 0

      # data[data == netCDF4.default_fillvals['f4']] = np.nan

      return data

   def _copy(self, data):
      data = data[:].astype(float)
      q = copy.deepcopy(data)
      # Remove missing values. Convert to -999 and then back to nan to avoid
      # warning messages when doing <, >, and == comparisons with nan.
      q[np.isnan(q)] = -999
      q[(q == -999) | (q < -1000000) | (q > 1e30)] = np.nan
      # q[q == netCDF4.default_fillvals["f4"]] = np.nan
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
      self.altitudes = np.zeros([1, 1])
      var_x = wxgen.variable.Variable("X")
      var_y = wxgen.variable.Variable("Y")
      var_z = wxgen.variable.Variable("Z")
      self.variables = [var_x, var_y, var_z]
      self._initial_state = {var_x: -10, var_y: -10, var_z: -25}
      self._std_initial_state = 10  # Standard deviation of initial condition error
      self.inittimes = np.ones(self.num) + wxgen.util.date_to_unixtime(20150101)
      self.timestep = 86400

      # Initialize
      self._data = dict()
      for var in self.variables:
         self._data[var] = np.zeros([self.length, self.X, self.Y, self.num], 'float32')
         self._data[var][0, 0, 0, :] = self._initial_state[var] + np.random.randn(self.num) * self._std_initial_state

      TT = int(1 / self._dt / 10)
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
      self.altitudes = np.nan * self.lats

      self._variance = 1
      self.timestep = 86400

   def _load(self, variable):
      values = np.random.randn(self.length, self.Y, self.X, self.num) * np.sqrt(self._variance)
      values = np.cumsum(values, axis=0)
      for i in range(self.length):
         # Ensure that the signal has a constant variance over time
         scale = 1./np.sqrt(1 + i)
         values[i, :, :, :] = values[i, :, :, :] * scale
      return values
