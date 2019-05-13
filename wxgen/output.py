import inspect
import netCDF4
import numpy as np
import sys
import os
import wxgen.util
import datetime


def get_all():
   """ Returns a list of all output classes """
   temp = inspect.getmembers(sys.modules[__name__], inspect.isclass)
   return temp


def get(name):
   """ Returns an output object of a class with the given name """
   outputs = get_all()
   m = None
   for mm in outputs:
      if(name == mm[0].lower()):
         m = mm[1]
   if m is None:
      wxgen.util.error("Cannot find output called '%s'" % name)
   return m


class Output(object):
   """
   A class for outputing trajectory information to file
   """
   def __init__(self, filename):
      self.filename = filename
      self.which_vars = None
      self.lat = None
      self.lon = None
      self.acc = None
      self.altitude = None
      self.write_indices = False
      self.command = None

   def write(self, trajectories, database, start_unixtime):
      """ Writes trajectories to file

      Arguments:
         trajectories (list): List of wxgen.trajectory
         database (wxgen.database): Database belonging to trajectories
         start_unixtime (int): Starting date for the output timeseries (unixtime)
      """
      raise NotImplementedError()


class Netcdf(Output):
   """
   Writes the trajectories to a netcdf file.
   """
   def write(self, trajectories, database, start_unixtime=wxgen.util.date_to_unixtime(20170101)):
      if len(trajectories) == 0:
         wxgen.util.error("No trajectories to write")
      file = netCDF4.Dataset(self.filename, 'w')

      file.createDimension("time")
      file.createDimension("ensemble_member", len(trajectories))
      use_single_gridpoint = self.lat is not None and self.lon is not None
      has_single_spatial_dim = database.X == 1
      xname = "longitude"
      yname = "latitude"

      if has_single_spatial_dim:
         file.createDimension('grid_point', database.Y)
         spatial_dims = ['grid_point']
      elif use_single_gridpoint:
         file.createDimension(yname, 1)
         file.createDimension(xname, 1)
         spatial_dims = [yname, xname]
      else:
         file.createDimension(yname, database.Y)
         file.createDimension(xname, database.X)
         spatial_dims = [yname, xname]

      # Time
      var_time = file.createVariable("time", "f8", ("time"))
      end_unixtime = start_unixtime + database.timestep * trajectories[0].length
      # TODO: This probably isn't right
      var_time[:] = np.arange(start_unixtime, end_unixtime, database.timestep)
      var_time.units = "seconds since 1970-01-01 00:00:00 +00:00"
      var_time.standard_name = "time"
      var_time.long_name = "time"

      # Forecast reference time
      var_frt = file.createVariable("forecast_reference_time", "f8")
      var_frt[:] = start_unixtime
      var_frt.units = "seconds since 1970-01-01 00:00:00 +00:00"
      var_frt.standard_name = "time"
      var_frt.long_name = "time"

      # Projection
      var_proj = file.createVariable("projection_regular_ll", "i4")
      var_proj.grid_mapping_name = "latitude_longitude"
      var_proj.earth_radius = 6367470
      var_proj.proj4 = "+proj=longlat +a=6367470 +e=0 +no_defs"

      if database.x is not None and database.y is not None:
         if has_single_spatial_dim:
            var_x = file.createVariable("x", "f8", ['grid_point'])
            var_y = file.createVariable("y", "f8", ['grid_point'])
            # var_z = file.createVariable("z", "f8", ['grid_point'])
         else:
            var_x = file.createVariable("x", "f8", xname)
            var_y = file.createVariable("y", "f8", yname)
            # var_z = file.createVariable("z", "f8", ('yname', 'xname'))
         var_x.units = database.x.units
         var_x.axis = database.x.axis
         var_x.standard_name = database.x.standard_name

         var_y.units = database.y.units
         var_y.axis = database.y.axis
         var_y.standard_name = database.y.standard_name

         # var_z.units = database.z.units
         # var_z.axis = database.z.axis
         # var_z.standard_name = database.z.standard_name

         var_x[:] = database.x[:]
         var_y[:] = database.y[:]
         # var_z[:] = database.z[:]

      if database.crs is not None:
         var_crs = file.createVariable("crs", "i4", [])
         var_crs.grid_mapping_name = database.crs.grid_mapping_name
         var_crs.proj4 = database.crs.proj4
         var_crs.epsg_code = database.crs.epsg_code

      # Latitude
      if has_single_spatial_dim:
         # Assume a lat/lon grid
         var_lat = file.createVariable("latitude", "f4", ('grid_point',))
         var_lat.units = "degrees_north"
         var_lat.standard_name = "latitude"
         if use_single_gridpoint:
            var_lat[:] = self.lat
         else:
            var_lat[:] = database.lats[:]

         # Longitude
         var_lon = file.createVariable("longitude", "f4", ('grid_point',))
         var_lon.units = "degrees_east"
         var_lon.standard_name = "longitude"
         if use_single_gridpoint:
            var_lon[:] = self.lon
         else:
            var_lon[:] = database.lons[:]

         # ATODO: Allow this to be called altitude, to conform to MET Norway standards
         var_altitude = file.createVariable("z", "f4", ('grid_point',))
         var_altitude.units = "m"
         var_altitude.standard_name = "height"
         if use_single_gridpoint:
            if self.altitude is not None:
               var_altitude[:] = self.altitude
         else:
            var_altitude[:] = database.altitudes[:]
      else:
         # Assume a lat/lon grid
         var_lat = file.createVariable("latitude", "f4", (yname))
         var_lat.units = "degrees_north"
         var_lat.standard_name = "latitude"
         if use_single_gridpoint:
            var_lat[:] = self.lat
         else:
            var_lat[:] = database.lats[:, 0]

         # Longitude
         var_lon = file.createVariable("longitude", "f4", (xname))
         var_lon.units = "degrees_east"
         var_lon.standard_name = "longitude"
         if use_single_gridpoint:
            var_lon[:] = self.lon
         else:
            var_lon[:] = database.lons[0, :]

         var_altitude = file.createVariable("z", "f4", (yname, xname))
         var_altitude.units = "m"
         var_altitude.standard_name = "height"
         if use_single_gridpoint:
            if self.altitude is not None:
               var_altitude[:] = self.altitude
         else:
            var_altitude[:] = database.altitudes[:]

      # Define forecast variables
      variables = database.variables
      vars = dict()
      for var in variables:
         if has_single_spatial_dim:
            vars[var.name] = file.createVariable(var.name, "f4", ("time", "ensemble_member", 'grid_point'))
         else:
            vars[var.name] = file.createVariable(var.name, "f4", ("time", "ensemble_member", yname, xname))
         if var.units is not None:
            vars[var.name].units = var.units
         if database.crs is not None:
            vars[var.name].grid_mapping = "crs"
            vars[var.name].coordinates = "y x z"

      # Write forecast variables
      if use_single_gridpoint:
         Xref, Yref = wxgen.util.get_i_j(database.lats, database.lons, self.lat, self.lon)
      for v in range(len(variables)):
         # Save variable after writing. Combine this loop with previous var loop.
         for m in range(len(trajectories)):
            values = database.extract_grid(trajectories[m], variables[v])
            # Insert a singleton dimension at dimension index 1
            values = np.expand_dims(values, 1)
            if use_single_gridpoint:
               vars[variables[v].name][:, m, :] = values[:, :, Xref, Yref]
            elif has_single_spatial_dim:
               vars[variables[v].name][:, m, :] = values[:, :, :, :]
            else:
               vars[variables[v].name][:, m, :, :] = values[:, :, :, :]
         if self.acc is not None and v in self.acc:
            vars[variables[v].name][:, ...] = np.cumsum(vars[variables[v].name], axis=0)

      if self.write_indices:
         var_segment_member = file.createVariable("segment_member", "i4", ("time", "ensemble_member"))
         var_segment_leadtime = file.createVariable("segment_lead_time", "i4", ("time", "ensemble_member"))
         var_segment_leadtime.units = "seconds"
         var_segment_time = file.createVariable("segment_time", "f8", ("time", "ensemble_member"))
         var_segment_time.units = "seconds since 1970-01-01 00:00:00 +00:00"
         var_segment_time.standard_name = "lead_time"
         var_segment_time.long_name = "lead_time"
         dt = database.timestep
         for m in range(0, len(trajectories)):
            trajectory = trajectories[m]
            var_segment_member[:, m] = trajectory.indices[:, 0]
            var_segment_leadtime[:, m] = database.leadtimes[trajectory.indices[:, 1]]
            var_segment_time[:, m] = database.inittimes[trajectory.indices[:, 0]]

      # Global attributes
      file.Conventions = "CF-1.0"
      history = "%s" % datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S +00:00')
      if self.command is None:
         history += ": wxgen"
      else:
         history += ": %s" % self.command
      file.history = history
      file.close()
