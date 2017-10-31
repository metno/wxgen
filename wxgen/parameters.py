import netCDF4
import numpy as np
import pyproj


class Parameters():
   """
   Class representing downscaling parameters. Implemented by a NetCDF file.

   Must contain
      - latitude (X, Y)
      - lonitude (X, Y)
      - proj

      One or more variables to represent the downscaling
   """
   def __init__(self, filename):
      self.filename = filename
      self.file = netCDF4.Dataset(self.filename, 'r')
      proj = pyproj.Proj(self.proj)
      if "latitude" not in self.file.variables or "longitude" not in self.file.variables:
         x, y = np.meshgrid(self.x, self.y)
         self._lons, self._lats = proj(x, y, inverse=True)
      else:
         self._lats = self.file.variables["latitude"][:]
         self._lons = self.file.variables["longitude"][:]
      self._field_cache = dict()

   @property
   def x(self):
      return self.file.variables["x"][:]

   @property
   def y(self):
      return self.file.variables["y"][:]

   @property
   def lats(self):
      return self._lats

   @property
   def lons(self):
      return self._lons

   @property
   def proj(self):
      """
      Returns:
         str: Proj4 string
      """
      return self.file.variables["projection_laea"].proj4

   def field(self, name):
      """
      Arguments:
         name (str): Name of field

      Returns:
         np.array: 3D array with dimensions X, Y, N
      """
      if name not in self._field_cache:
         values = self.file.variables[name][:].astype(float)
         values = np.ma.filled(values, fill_value=np.nan)
         self._field_cache[name] = values
      return self._field_cache[name]
