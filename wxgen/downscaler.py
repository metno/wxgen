import inspect
import numpy as np
import sys
import wxgen.util
import scipy.spatial


def get_all():
   """ Returns a list of all metric classes """
   temp = inspect.getmembers(sys.modules[__name__], inspect.isclass)
   return temp


def get(name):
   """ Returns a downscaler object of a class with the given name """
   downscalers = get_all()
   m = None
   for mm in downscalers:
      if(name == mm[0].lower()):
         m = mm[1]()
   if m is None:
      wxgen.util.error("Cannot find downscaler called '%s'" % name)
   return m


class Downscaler(object):
   """
   Downscale large-scale trajectories
   """
   def __init__(self):
      pass

   def generate(self, values):
      """

      Arguments:
         values (np.array): A large-scale trajectory

      Returns:
         np.array: A 4D array (T, X, Y, V)
      """
      raise NotImplementedError()


class Qq(object):
   """
   Quantile mapping between a fine-scale model and a coarse-scale model
   """
   def __init__(self):
      self.cache_nn = None

   def generate(self, raw, lats, lons, parameters):
      # Where should the nearest neighbour stuff be done?
      X = parameters.lons.shape[1]
      Y = parameters.lons.shape[0]
      v0 = parameters.field("coarse_scale")
      v1 = parameters.field("fine_scale")
      values = np.zeros([Y, X])
      lats_d = parameters.lats
      lons_d = parameters.lons

      # Downscale the raw forecast
      coords = np.zeros([len(lats.flatten()), 2])
      coords[:, 0] = lats.flatten()
      coords[:, 1] = lons.flatten()
      if self.cache_nn is None:
         tree = scipy.spatial.KDTree(coords)
         self.cache_nn = np.zeros([Y, X], 'int')
         for y in range(Y):
            for x in range(X):
               dist, index = tree.query([lats_d[y, x], lons_d[y, x]])
               self.cache_nn[y, x] = index
      for y in range(Y):
         for x in range(X):
            values[y, x] = np.interp(raw.flatten()[self.cache_nn[y, x]], v0[:, y, x], v1[:, y, x])
      # values = np.random.rand(Y, X)
      return values


class NearestNeighbour(object):
   """
   Nearest neighbour
   """
   def __init__(self):
      self.cache_nn = None

   def generate(self, raw, lats, lons, parameters):
      # Where should the nearest neighbour stuff be done?
      X = parameters.lons.shape[1]
      Y = parameters.lons.shape[0]
      v0 = parameters.field("coarse_scale")
      v1 = parameters.field("fine_scale")
      values = np.zeros([Y, X])
      lats_d = parameters.lats
      lons_d = parameters.lons

      # Downscale the raw forecast
      coords = np.zeros([len(lats.flatten()), 2])
      coords[:, 0] = lats.flatten()
      coords[:, 1] = lons.flatten()
      if self.cache_nn is None:
         tree = scipy.spatial.KDTree(coords)
         self.cache_nn = np.zeros([Y, X], 'int')
         for y in range(Y):
            for x in range(X):
               dist, index = tree.query([lats_d[y, x], lons_d[y, x]])
               self.cache_nn[y, x] = index
      for y in range(Y):
         for x in range(X):
            values[y, x] = raw.flatten()[self.cache_nn[y, x]]
      return values
