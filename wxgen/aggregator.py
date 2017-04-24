import numpy as np
import inspect
import sys
import wxgen.util


def get_all():
   """ Returns a list of all aggregator classes """
   temp = inspect.getmembers(sys.modules[__name__], inspect.isclass)
   return temp


def get(name):
   """
   Returns an instance of an object with the given class name

   Arguments:
      name (str): The name of the class. Use a number between 0 and 1 to get the corresponding
         quantile aggregator.
   """
   aggregators = get_all()
   a = None
   for aggregator in aggregators:
      if(name == aggregator[0].lower()):
         a = aggregator[1]()
   if a is None and wxgen.util.is_number(name):
      a = Quantile(float(name))

   return a


class Aggregator(object):
   """
   Base class for aggregating an array (i.e. computing a single value from an
   array)

   Usage:
   mean = wxgen.aggregator.Mean()
   mean(np.array([1,2,3]))
   """
   def __hash__(self):
      # TODO
      return 1

   def __eq__(self, other):
      return self.__class__ == other.__class__

   @classmethod
   def name(cls):
      """ Returns a string representing the name of the aggregator """
      return cls.__name__.lower()

   def __call__(self, array, axis=None):
      """ Implement this function to return the scalar value """
      raise NotImplementedError()


class Mean(Aggregator):
   def __call__(self, array, axis=None):
      return np.mean(array, axis=axis)


class Median(Aggregator):
   def __call__(self, array, axis=None):
      return np.median(array, axis=axis)


class Min(Aggregator):
   def __call__(self, array, axis=None):
      return np.min(array, axis=axis)


class Max(Aggregator):
   def __call__(self, array, axis=None):
      return np.max(array, axis=axis)


class Std(Aggregator):
   def __call__(self, array, axis=None):
      N = len(array)
      return np.std(array, axis=axis, ddof=1)


class Variance(Aggregator):
   def __call__(self, array, axis=None):
      return np.var(array, axis=axis, ddof=1)


class Iqr(Aggregator):
   def __call__(self, array, axis=None):
      return np.percentile(array, 75, axis=axis) - np.percentile(array, 25, axis=axis)


class Range(Aggregator):
   def __call__(self, array, axis=None):
      return wxgen.util.nprange(array=axis)


class Sum(Aggregator):
   def __call__(self, array, axis=None):
      return np.sum(array, axis=axis)


class Consecutive(Aggregator):
   def __init__(self):
      self.threshold = 1

   def __call__(self, array, axis=None):
      acc = np.cumsum(array == self.threshold, axis=axis)

      if axis is None:
         last_acc = 0
         for i in range(1, len(acc)):
            if array[i] != self.threshold:
               last_acc = acc[i]
               acc[i] = 0
            else:
               acc[i] = acc[i] - last_acc
      elif axis == 0:
         shape = np.array(array.shape)
         shape = np.delete(shape, 0)
         max = np.zeros(shape)
         acc = np.zeros(shape)
         for i in range(array.shape[0]):
            acc = (acc + array[i, ...])*array[i, ...]
            # Use np.maximum not np.max, as this will compare two separate arrays
            max = np.maximum(max, acc)
         return max
      elif axis == 1:
         shape = np.array(array.shape)
         shape = np.delete(shape, 1)
         max = np.zeros(shape)
         acc = np.zeros(shape)
         for i in range(array.shape[1]):
            acc = (acc + array[:, i, ...])*array[:, i, ...]
            # Use np.maximum not np.max, as this will compare two separate arrays
            max = np.maximum(max, acc)
         return max
      else:
         raise NotImplementedError()

      if 1:
         pass
      elif axis == 0:
         if len(acc.shape) == 1:
            acc = np.expand_dims(acc, 1)
            array = np.expand_dims(array, 1)
         for j in range(acc.shape[1]):
            last_acc = 0
            for i in range(1, acc.shape[0]):
               if array[i, j] != self.threshold:
                  last_acc = acc[i, j]
                  acc[i, j] = 0
               else:
                  acc[i, j] = acc[i, j] - last_acc
      elif axis == 1:
         if len(acc.shape) == 1:
            acc = np.expand_dims(acc, 1)
            array = np.expand_dims(array, 1)
         for i in range(acc.shape[0]):
            last_acc = 0
            for j in range(1, acc.shape[1]):
               if array[i, j] != self.threshold:
                  last_acc = acc[i, j]
                  acc[i, j] = 0
               else:
                  acc[i, j] = acc[i, j] - last_acc

      return np.max(acc, axis=axis)


class Quantile(Aggregator):
   """ Returns a specified quantile from an array """
   def __init__(self, quantile):
      self.quantile = quantile
      if self.quantile < 0 or self.quantile > 1:
         wxgen.util.error("Quantile must be between 0 and 1")

   def __call__(self, array):
      return np.percentile(array, self.quantile*100)

   def __eq__(self, other):
      return (self.__class__ == other.__class__) and (self.quantile == other.quantile)
