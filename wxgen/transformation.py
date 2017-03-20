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
      wxgen.util.error("Cannot find transformation called '%s'" % name)
   return m


class Transformation(object):
   def transform(self, array):
      """ Transforms values in an array

      Arguments:
         array (np.array): Array of values

      Returns:
         np.array: Array of transformed values
      """
      raise NotImplementedError()


class Nothing(Transformation):
   """ No transformation """
   def transform(self, array):
      return array


class FrostDay(Transformation):
   """ 1 if temperature is freezing """
   threshold = 273.15

   def transform(self, array):
      return array < self.threshold


class SummerDay(Transformation):
   """ 1 if temperature is above 10 C """
   threshold = 273.15 + 10

   def transform(self, array):
      return array > self.threshold


class DryDay(Transformation):
   """ 1 if precip is less than 1 """
   threshold = 1

   def transform(self, array):
      return array < self.threshold
