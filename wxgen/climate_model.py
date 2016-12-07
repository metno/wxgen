import numpy as np
import datetime

def get(name):
   if name == "30":
      return Bin(30)

class ClimateModel(object):
   def get(self, unixtimes):
      """ Returns a representation of the state for a given date
      
      Arguments:
      unixtimes       A numpy array of unix times
      """
      raise NotImplemented


class Bin(ClimateModel):
   def __init__(self, num_days):
      self._num_days = num_days

   def get(self, unixtimes):
      day_of_year = np.array([int(datetime.datetime.fromtimestamp(unixtime).strftime('%j')) for
         unixtime in unixtimes])
      #day_of_year = np.array([datetime.datetime.utcfromtimestamp(int(unixtime)).timetuple().tm_yday for
      #      unixtime in unixtimes])
      return day_of_year / self._num_days


class Index(ClimateModel):
   def __init__(self, file, num_days=365):
      self._file = file
      self._num_days = num_days

   def get(self, date):
      return np.nan
