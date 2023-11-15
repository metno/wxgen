import copy
import os
import time as timing
from typing import TypeAlias

import netCDF4
import numpy as np
import pywt
import xarray as xr

import wxgen.climate_model
import wxgen.config
from wxgen.trajectory import Trajectory
import wxgen.util
import wxgen.variable

import logging


class Database:
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
       x (np.array): 1D array of x-coordinates
       y (np.array): 1D grid of y-coordinates
       crs (np.array): Projection variable
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
        self.join_config_fn = None
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
        self.deacc = None
        self.x = None
        self.y = None
        # self.z = None
        self.crs = None

        self._join_config_cache = None

        self.logger = logging.getLogger(self.__class__.__name__)

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
            self.logger.debug("Cache miss variable '%s'", variable.name)

            if len(self._data_cache) > 0:
                # Check if we need to remove data from cache

                """
                Calculate the size of the data in the cache. Explicitly convert shape to int64,
                since on Windows64 this defaults to int32 and can causes overflow in np.product.
                """
                akey = list(self._data_cache.keys())[0]
                bytes_per_value = 4
                shape = np.array(self._data_cache[akey].shape, 'int64')
                size_per_key = np.product(shape) * bytes_per_value

                next_size = float(len(self._data_cache) + 1) * size_per_key
                next_size_gb = next_size / 1e9
                if self.mem is not None and next_size_gb > self.mem:
                    # remove from cache
                    I = np.random.randint(len(self._data_cache))
                    rmkey = list(self._data_cache.keys())[I]
                    self._data_cache.pop(rmkey)
                    self.logger.warning("Cache full (%2.1fGB): Removing member '%s' from cache", next_size_gb, rmkey.name)

            # Get data from subclass
            data = self._load(variable)
            if self.deacc is not None and variable.name in self.deacc:
                data[1:, ...] = np.diff(data, axis=0)
                I = np.isnan(data[0, ...])
                data[0, ...] = 0
                # Preserve nans
                data[0, ...][I] = np.nan
                self.logger.debug("Deaccumulating variable: %s", variable.name)
            self._data_cache[variable] = data
            e = timing.time()
            # print "Timing: %g" % (e - t)
        else:
            self.logger.debug("Cache hit '%s'", variable.name)

        return self._data_cache[variable]

    def _load(self, variable):
        raise NotImplementedError()

    def get(self, i: int, i_lead_time_start: int = 0) -> Trajectory:
        """ Get the i'th trajectory in the database.

        Args:
            i: segment index
            i_lead_time_start: first index of the lead_times to consider
        """
        indices = np.zeros([self.length-i_lead_time_start, 2], int)
        indices[:, 0] = i
        indices[:, 1] = np.arange(i_lead_time_start, self.length)
        assert(np.sum(np.isnan(indices)) == 0)
        return Trajectory(indices)

    @property
    def remove_first_timestep(self):
        return self.deacc is not None

    def get_truth(self, start_date=None, end_date=None, start_time_of_day=None):
        """
        Concatenate the initialization state of all trajectories

        Arguments:
           start_date (int): Starting date of scenario (YYYYMMDD)
           end_date (int): End date of scenario (YYYYMMDD)
           start_time_of_day (int): What time of day (in seconds) should the scenario start and end at?

        Returns:
           trajectory.Trajectory: Truth trajectory

        Note that if (for example) start_date=20180101, end_date=20180103, start_time_of_day=0, then the
        tracjectory ends at 00Z on 20180103, i.e. it doesn't include the whole end_date.

        """
        if start_time_of_day is not None and start_time_of_day not in self.leadtimes:
            RuntimeError("A start time of %d h is not possible when the timestep is %d h" % (start_time_of_day // 3600, self.timestep // 3600))

        """ Compute the unixtimes to find truth for """
        include_extra_time_step_at_end = True
        start = np.min(self.inittimes) if start_date is None else wxgen.util.date_to_unixtime(start_date)
        end = np.max(self.inittimes) if end_date is None else wxgen.util.date_to_unixtime(end_date)
        if include_extra_time_step_at_end:
            times = np.arange(start, end + self.timestep, self.timestep)
        else:
            times = np.arange(start, end, self.timestep)

        indices = -1*np.ones([len(times), 2], int)
        self.logger.debug("Start: %d End: %d", start, end)
        for i in range(len(times)):
            time = times[i]

            """
            When we have accumulatd variables, we allow the first timestep to be used at the
            beginning of the scenario. However, later in the scenario we cannot use it.
            """
            if i == 0 or not self.remove_first_timestep:
                I = np.where(time >= self.inittimes)[0]
            else:
                I = np.where(time > self.inittimes)[0]
            if len(I) == 0:
                RuntimeError("There are no inittimes available before %d. The earliest is %d." % (wxgen.util.unixtime_to_date(time), wxgen.util.unixtime_to_date(np.min(self.inittimes))))

            """ Find the inittime and leadtime that is closest to the desired date """
            curr_times = self.inittimes[I]
            Ibest = np.argmin(np.abs(curr_times - time))
            inittime = self.inittimes[I][Ibest]
            lt = int((time - inittime)/self.timestep)
            if lt < self.length:
                indices[i, 0] = np.where(self.inittimes == inittime)[0][0]
                indices[i, 1] = lt
            else:
                assert(i > 0)
                dI = int(86400 / self.timestep)
                if i - dI < 0:
                    RuntimeError("It is not possible to create a truth scenario because there are missing days and the segments are too short to cover the holes. Processing stopped at %dT%d." % (wxgen.util.unixtime_to_date(time), wxgen.util.unixtime_to_hour(time)))
                self.logger.warning("Did not find an index for %d = %dT%02dZ. Using the previous day's state.", time, wxgen.util.unixtime_to_date(time), wxgen.util.unixtime_to_hour(time))
                indices[i, 0] = indices[i-dI, 0]
                indices[i, 1] = indices[i-dI, 1]

        if start_time_of_day is not None:
            """ Crop the timeseries such that the starting and ending hour of the day is as desired """
            Iindices = np.where(self.leadtimes[indices[:, 1]] % 86400 == start_time_of_day)[0][0]
            indices = indices[Iindices:, :]

            """ Crop away at the end """
            Iindices = np.where(self.leadtimes[indices[:, 1]] % 86400 == start_time_of_day)[0][-1]
            if include_extra_time_step_at_end:
                indices = indices[:(Iindices+1), :]
            else:
                indices = indices[:Iindices, :]

            # Crop last few hours
            # TODO: Deal with the fact that the cropping makes it so that the trajectory is too short
            if start_time_of_day != 0:
                self.logger.warning("The scenario will likely be one day too short, since the requested start time is not 00Z. A fix for this has not been implemented yet.")

        return Trajectory(indices)

    def extract(self, trajectory: Trajectory) -> np.ndarray:
        """
        Extract a trajectory of large-scale aggregated values from the database

        Arguments:
           trajectory (Trajectory): Trajectory to extract

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

    def extract_grid(self, trajectory: Trajectory, variable):
        """
        Extract a trajectory of large-scale values from the database

        Arguments:
           trajectory (Trajectory): Trajectory to extract
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

    def extract_matching(self, trajectory: Trajectory):
        """
        Extract a trajectory of values used to match states from the database

        Arguments:
           trajectory: Trajectory to extract

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
    def join_config(self) -> wxgen.config.Config:
        if self._join_config_cache is None:
            self._join_config_cache = wxgen.config.Config(self.join_config_fn)
        return self._join_config_cache

    @property
    def _data_matching(self) -> np.ndarray:
        """Data used to match states
        
        This may be different than data_agg, since it can include wavelet information.

        Returns:
           Array with dimensions [segment_lead_time, variable_at_point, idx_segment], where `variable_at_point` are
           the variable values at a the matching points `p`.
        """
        if self._data_matching_cache is None:
            if self.join_config_fn is not None:
                # Preload all required variables
                self._data_matching_cache = np.zeros([self.length, len(self.join_config.points), self.num], 'float32')

                count = 0
                for variable_name, points_df in self.join_config.group_by_variable().items():
                    variable = self.get_variable_by_name(variable_name)
                    if variable == None:
                        RuntimeError("Cannot use variable %s in join config, since this is not in the database" % variable_name)
                    temp = self.load(variable)
                    for _, point in points_df.iterrows():
                        # Find nearest neighbour
                        dist = (self.lats - point['lat']) ** 2 + (self.lons - point['lon']) ** 2
                        indices = np.unravel_index(dist.argmin(), dist.shape)
                        self._data_matching_cache[:, count, :] = temp[:, indices[0], indices[1], :]
                        count += 1

            elif self.spatial_decomposition == 0:
                self._data_matching_cache = self._data_agg
            elif self.spatial_decomposition == 'all':
                N = int(self.X * self.Y)
                size = [self.length, self.V*N, self.num]
                self.logger.debug("Allocating %g GB for matching values", (np.product(size) * 4.0 / 1e9))

                self._data_matching_cache = np.zeros(size, 'float32')
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
                self.logger.debug("Computing wavelets")
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
                self.logger.debug("Wavelet time: %f", (e - s))

        return self._data_matching_cache

    @property
    def leadtimes(self):
        return np.arange(0, (self.length+1) * self.timestep, step=self.timestep)

    """
        Returns the variable object with the given variable name
    """
    def get_variable_by_name(self, name):
        for variable in self.variables:
            if variable.name == name:
                return variable
        return None


SegmentIndices: TypeAlias = np.ndarray
"""Array with index of sements in the database.

Can e.g. be used to collect the indices of all segments that have valid values (non-nan for all variables)
"""


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
            RuntimeError("File '%s' does not exist" % filename)

        self.name = filename[filename.rfind('/') + 1:]
        self.filename = filename

        self.mem = mem
        self._file = netCDF4.Dataset(filename)

        # Set dimensions
        var_names = [name for name in self._file.variables if name not in ["lat", "lon", "latitude",
           "longitude", "altitude", "x", "y", "z", "ensemble_member", "lead_time", "dummy",
           "forecast_is_complete", "crs",
           "longitude_latitude", "time", "projection_regular_ll", "segment_lead_time", "segment_member", "segment_time"]]
        if vars is None:
            vars = range(len(var_names))

        """ Check if all variables are available """
        for var in vars:
            if wxgen.util.is_number(var):
                if var >= len(var_names):
                    RuntimeError("Input has only %d variables. Variable %d is outside range." % (len(var_names), max(vars)))
            else:
                if var not in var_names:
                    RuntimeError("Input does not have variable '%s'" % var)

        if "lead_time" in self._file.dimensions and "time" in self._file.dimensions:
            lead_time_dim = "lead_time"
            time_dim = "time"
        else:
            RuntimeError("Cannot read file. Missing 'time' and/or 'lead_time' dimensions")
        self.variables = list()
        for i, var in enumerate(vars):
            if wxgen.util.is_number(var):
                var_name = var_names[var]
            else:
                var_name = var
            units = None
            label = None
            dims = self._file.variables[var_name].dimensions
            if "time" in dims:
                if hasattr(self._file.variables[var_name], "units"):
                    units = self._file.variables[var_name].units
                if hasattr(self._file.variables[var_name], "standard_name"):
                    label = self._file.variables[var_name].standard_name
                var = wxgen.variable.Variable(var_name, units, label)
                if var in self.variables:
                    self.logger.warning("Variable '%s' is specified multiple times. Only using it once.", var.name)
                else:
                    self.variables.append(var)
                    self.logger.debug("Using variable '%s'", var_name)

        # Load data
        lead_time_var = self._file.variables[lead_time_dim]

        # TODO: SEEMS LIKE self.length is set wrong, if if should be # days! - why is it not reset?
        # Assume dataset has a daily timestep, except if it is possible to deduce otherwise
        self.timestep = 86400
        leadtimes = wxgen.util.clean(self._file.variables[lead_time_dim][:])
        self.length = len(leadtimes)
        if hasattr(lead_time_var, "units"):
            if len(lead_time_var.units) >= 7 and lead_time_var.units[0:7] == "seconds":
                pass
            elif lead_time_var.units == "days":
                leadtimes *= 86400
            elif lead_time_var.units == "hours":
                leadtimes *= 3600
            else:
                wxgen.error("Cannot parse time units")
        self.first_lead_time = leadtimes[0]

        if len(np.unique(np.diff(leadtimes))) != 1:
            RuntimeError("Cannot handle databases with leadtimes that are not spaced evenly: %s" %
                  (','.join(["%g" % leadtime for leadtime in leadtimes])))

        self.timestep = leadtimes[1] - leadtimes[0]
        # # TODO: at least approx -- should find out on whether to use ceil or floor...
        # self.length = int(len(leadtimes) / (86400 / self.timestep)) 

        if np.isnan(self.timestep):
            RuntimeError("Cannot determine timestep from database")

        self.has_single_spatial_dim = False
        times = self._file.variables["time"][:]
        self.ens = len(self._file.dimensions["ensemble_member"])

        """ Find for which time indices there is valid data """
        if "forecast_is_complete" in self._file.variables:
            self.logger.debug("Using forecast_is_complete to remove missing times")
            self._Itimes = np.where((np.isnan(times) == 0) & (self._file.variables["forecast_is_complete"][:] == 1))[0]
            if len(self._Itimes) != len(times):
                self.logger.debug("Removing %d times due to missing values", (len(times) - len(self._Itimes)))
        else:
            self._Itimes = np.where(np.isnan(times) == 0)[0]
        times = times[self._Itimes]
        self.inittimes = np.repeat(times, self.ens)
        self.num = len(self._Itimes) * self.ens

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
        self.logger.debug("Has single spatial dimension? %d", self.has_single_spatial_dim)

        """
        Read projection information
        """
        if "x" in self._file.variables:
            self.x = self._file.variables['x']
        if "y" in self._file.variables:
            self.y = self._file.variables['y']
        # if "z" in self._file.variables:
        #    self.z = self._file.variables['z']
        if "crs" in self._file.variables:
            self.crs = self._file.variables['crs']

        """
        Read lat/lon/elev variables
        """
        if self.has_single_spatial_dim:
            if "latitude" in self._file.variables:
                self.lats = self._copy(self._file.variables["latitude"])
                self.lons = self._copy(self._file.variables["longitude"])
            elif "lat" in self._file.variables:
                self.lats = self._copy(self._file.variables["lat"])
                self.lons = self._copy(self._file.variables["lon"])
            else:
                if self.x is not None and self.y is not None:
                    self.logger.debug("Diagnosing lat/lon from projection information")
                    self.lats, self.lons = wxgen.util.get_latlon_from_proj(self.crs.proj4, self.x, self.y)
                    #proj = pyproj.Proj(self.crs.proj4)
                    #self.lons, self.lats = proj(self.x[:], self.y[:], inverse=True)
                else:
                    self.logger.warning("Could not determine lat/lon values")
                    self.lons = np.nan * np.zeros(len(self.x))
                    self.lats = np.nan * np.zeros(len(self.x))

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
                self.logger.debug("Meshing latitudes and longitudes")
                self.lons, self.lats = np.meshgrid(self.lons, self.lats)
        """
        Read altitude information
        """
        if "altitude" in self._file.variables:
            self.altitudes = self._copy(self._file.variables["altitude"])
        elif "z" in self._file.variables:
            self.altitudes = self._copy(self._file.variables["z"])
        elif "surface_geopotential" in self._file.variables:
            self.altitudes = np.squeeze(self._copy(self._file.variables["surface_geopotential"]/9.81))
        else:
            self.altitudes = np.nan * self.lats

        if self.altitudes.shape != self.lats.shape:
            # Try to remove singleton dimensions
            self.altitudes = np.squeeze(self.altitudes)

        if self.has_single_spatial_dim:
            self.altitudes = np.expand_dims(self.altitudes, 1)
            if len(self.altitudes.shape) == 1:
                self.altitudes = np.expand_dims(self.altitudes, 1)

        if self.altitudes.shape != self.lons.shape:
            RuntimeError("Altitude dimensions do not match those of lat/lon")

        self._src_time = None
        self._src_ensemble_member = None      


    def _load(self, variable) -> np.ndarray:
        """Reads from netcdf

        It processes data from an input form 
           [forecast issue time, lead_time, Y, X, ensemble member], or
           [forecast issue time, lead_time, grid_point, ensemble member]
        to 
           [lead_time, Y, X, segment_member],

        where the segment_member is extracted from the product of forecast issue time and the ensemble_member.
        If there is a only a single spatial dimension (grid-point), the data is mapped to two dimensions, X, Y, but
        one of them has length=1. 
        
        Returns:
            # TODO: [self.length, self.Y, self.X, self.num] or [self.length, self.X, self.Y, self.num]
            Data with dimensions [lead_time, Y, X, segment_member]
        """
        if variable.name not in self._file.variables:
            RuntimeError("Variable '%s' does not exist in file '%s'" % (variable.name, self.name))
        self.logger.debug("Allocating %g GB for '%s'", np.product(self._file.variables[variable.name].shape) * 4.0 / 1e9, variable.name)

        # Loading whole data-set and perfrom slicing afterwards is much faster. 
        # Additionally, h5netcdf seemed to be ~2x faster then default.
        with xr.open_dataset(self.filename, engine='h5netcdf') as data_file:
            data = data_file[variable.name].load()
            # src_time = data_file.time.load()
            # ensemble_member = data_file.ensemble_member.load()

        # remove where not 'forecast_is_complete'
        data = data[self._Itimes, ...]
        # src_time = src_time[self._Itimes]

        # same action, if forecast_is_complete did not exists
        if np.isnan(data).any():
            arr_data_is_finite = np.isfinite(data).all(dim=["lead_time", "grid_point", "ensemble_member"])
            idx_time_is_finite = np.nonzero(arr_data_is_finite.values)[0]
            self.logger.warning(f"Removing {len(src_time) - len(idx_time_is_finite)}  member due to nan values.")
            # src_time = src_time[idx_time_is_finite]
            data = data.isel(time=idx_time_is_finite)

        data = data.stack(n_segment=["time", "ensemble_member"])
        self._src_time = data.time.values
        self._src_ensemble_member = data.ensemble_member.values

        if self.has_single_spatial_dim:
            data = data.expand_dims(dim={"X": 1}, axis=2)
            data = data.values.copy() # copy, since expand_dims is a view only
        else:
            data = data.values

        # Hotfix
        if variable.name == "precipitation_amount":
            data[data < 0] = 0

        # data[data == netCDF4.default_fillvals['f4']] = np.nan

        return data

    @property
    def src_time(self) -> None | np.ndarray["datetime64[ns]"]:
        """'time' variable of each segment_member in the database"""
        return self._src_time
        
    @property
    def src_ensemble_member(self) -> None | np.ndarray[int]:
        """'ensemble_member' variable of each segment_member in the database"""
        return self._src_ensemble_member

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
