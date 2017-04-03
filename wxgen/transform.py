import inspect
import numpy as np
import sys
import wxgen.util


def get_all():
   """ Returns a list of all classes """
   temp = inspect.getmembers(sys.modules[__name__], inspect.isclass)
   return temp


def get(name):
   """ Returns an object of a class with the given name """
   classes = get_all()
   m = None
   for mm in classes:
      if(name == mm[0].lower()):
         m = mm[1]()
   if m is None:
      wxgen.util.error("Cannot find transform called '%s'" % name)
   return m


class Transform(object):
   def __call__(self, array):
      """ Transforms values in an array

      Arguments:
         array (np.array): Array of values

      Returns:
         np.array: Array of transformed values
      """
      raise NotImplementedError()


class Nothing(Transform):
   """ No transform """
   def __call__(self, array):
      return array


class FrostDay(Transform):
   """ 1 if temperature is freezing """
   threshold = 273.15

   def __call__(self, array):
      return array < self.threshold


class SummerDay(Transform):
   """ 1 if temperature is above 10 C """
   threshold = 273.15 + 10

   def __call__(self, array):
      return array > self.threshold


class DryDay(Transform):
   """ 1 if precip is less than 1 """
   threshold = 1

   def __call__(self, array):
      return array < self.threshold


class WetDay(Transform):
   """ 1 if precip is greater than 1 """
   threshold = 1

   def __call__(self, array):
      return array >= self.threshold
