import numpy as np
import datetime
import wxgen.util

def get(name):
   if name == "30":
      return Bin(30)

def day_of_year(unixtimes):
   return np.array([int(datetime.datetime.fromtimestamp(unixtime).strftime('%j')) for unixtime in unixtimes])

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
      #day_of_year = np.array([datetime.datetime.utcfromtimestamp(int(unixtime)).timetuple().tm_yday for
      #      unixtime in unixtimes])
      day = day_of_year(unixtimes)
      return day / self._num_days


class Index(ClimateModel):
   def __init__(self, filename, num_days=365):
      self._filename = filename
      self._num_days = num_days
      self._edges = np.array([-100,-1, 1, 100])
      fid = open(self._filename, 'r')
      self._index = dict()
      prev = 0
      for line in fid:
         words = line.strip().split(' ')
         words = [word for word in words if word != ""]
         year = int(words[0])
         month = int(words[1])
         day = int(words[2])
         date = int("%04d%02d%02d" % (year, month, day))
         unixtime = wxgen.util.date_to_unixtime(date)
         if words[3] == "NA":
            value = prev
         else:
            value = float(words[3])
            prev = value
         I = np.where(value > self._edges)[0]
         self._index[unixtime] = I[-1]

   def get(self, unixtimes):
      day = day_of_year(unixtimes)/self._num_days
      index = np.nan*np.zeros(len(unixtimes))
      for i in range(0, len(unixtimes)):
         unixtime = unixtimes[i]
         if unixtime in self._index:
            index[i] = self._index[unixtime]
         else:
            print "Missing: %d" % unixtime
      print index
      return day + index * 100
