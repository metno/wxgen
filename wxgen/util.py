import numpy as np
import calendar
import datetime
import sys


DEBUG = False
COLORS = {"red": 31, "yellow": 33, "green": 32}


def random_weighted(weights):
   """
   Randomly selects an index into an array based on weights

   weights     A numpy array of weights
   """
   acc = np.cumsum(weights) / np.sum(weights)
   r = np.random.rand(1)
   temp = np.where(acc < r)[0]
   if len(temp) == 0:
      I = 0
   else:
      I = temp[-1]
   return I


def error(message):
   print "\033[1;31mError: " + message + "\033[0m"
   sys.exit(1)


def warning(message):
   print "\033[1;33mWarning: " + message + "\033[0m"


def debug(message, color="green"):
   col = COLORS[color]
   if DEBUG:
      print "\033[1;%imDebug: " % (col) + message + "\033[0m"


def parse_dates(dates):
   return [int(date) for date in parse_numbers(dates, True)]


def parse_ints(ints):
   return [int(num) for num in parse_numbers(ints, True)]


def parse_numbers(numbers, isDate=False):
   """
   Convert a string into an array of numbers. allowable formats:
      num
      num1,num2,num3
      start:end
      start:step:end
   """
   if numbers is None:
      return None

   # Check if valid string
   if(any(char not in set('-01234567890.:,') for char in numbers)):
      error("Could not translate '" + numbers + "' into numbers")

   values = list()
   commaLists = numbers.split(',')
   for commaList in commaLists:
      colonList = commaList.split(':')
      if(len(colonList) == 1):
         value = float(colonList[0])
         if int(value) == value:
            values.append(int(value))
         else:
            values.append(value)
      elif(len(colonList) <= 3):
         start = float(colonList[0])
         step = 1
         if(len(colonList) == 3):
            step = float(colonList[1])
         end = float(colonList[-1])
         if(isDate):
            date = min(start, end)
            curr = list()
            while date <= max(start, end):
               curr.append(date)
               date = get_date(date, step)
            values = values + list(curr)
         else:
            if int(start) == start and int(end) == end and int(step) == step:
               # Generate integer list
               values = values + list(range(int(start), int(end)+1, int(step)))
            else:
               # arange does not include the end point:
               stepSign = step / abs(step)
               values = values + list(np.arange(start, end + stepSign*0.0001, step))
      else:
         error("Could not translate '" + numbers + "' into numbers")
      if(isDate):
         for i in range(0, len(values)):
            values[i] = int(values[i])
   return values


def resize(vec, size):
   """
   Resizes a vector such that it has the right size. This is done by repeating the vector
   in each dimension until the required size is reached. Note an error is thrown if 'size'
   is not a multiple of the size of vec.

   vec      A 1D or 2D numpy array
   size     A list of dimension sizes (e.g. [2,3])
   """
   if not isinstance(vec, (np.ndarray)):
      vec_resized = vec * np.ones(size)
   elif vec.shape[0] == size[0] and len(vec.shape) == 1:
      vec_resized = np.reshape(np.repeat(vec, size[1]), size)
   elif vec.shape[0] == 1 and len(vec.shape) == 1:
      vec_resized = vec*np.ones(size)
   else:
      # Check that the output dims are multiples of input dims
      assert(size[0] % vec.shape[0] == 0)
      assert(size[1] % vec.shape[1] == 0)
      vec_resized = np.tile(vec, (size[0] / vec.shape[0], size[1] / vec.shape[1]))
   return vec_resized


def date_to_unixtime(date):
   year = date / 10000
   month = date / 100 % 100
   day = date % 100
   ut = calendar.timegm(datetime.datetime(year, month, day).timetuple())
   return ut


def unixtime_to_date(unixtime):
   dt = datetime.datetime.utcfromtimestamp(int(unixtime))
   date = dt.year * 10000 + dt.month * 100 + dt.day
   return date


def correlation(ar1, ar2, axis):
   """
   Computes the correlation between two multi-dimensional arrays
   """
   ar1mean = np.resize(np.mean(ar1, axis=axis), ar1.shape)
   ar2mean = np.resize(np.mean(ar2, axis=axis), ar2.shape)

   cov = np.mean((ar1 - ar1mean) * (ar2 - ar2mean), axis=axis)
   var1 = np.var(ar1, axis=axis)
   var2 = np.var(ar2, axis=axis)

   return cov / np.sqrt(var1) / np.sqrt(var2)


def get_date(date, diff):
   """ Date calculation: Adds 'diff' to 'date'

   Arguments:
      date (int): An integer of the form YYYYMMDD
      diff (int): Number of days to add to date

   Returns:
      int: A new date in the form YYYYMMDD
   """
   year = int(date / 10000)
   month = int(date / 100 % 100)
   day = int(date % 100)
   date2 = datetime.datetime(year, month, day, 0) + datetime.timedelta(diff)
   return int(date2.strftime('%Y%m%d'))


def day_of_year(unixtimes):
   """
   Arguments:
      unixtime (int): Number of seconds since 1970-01-01

   Returns:
      int: Day of year
   """
   ar = np.zeros([len(unixtimes), 1], int)
   ar[:, 0] = [int(datetime.datetime.fromtimestamp(unixtime).strftime('%j')) for unixtime in unixtimes]
   return ar


def get_i_j(lats, lons, lat, lon):
   """
   Finds the nearest neighbour in a lat lon grid. If the point is outside the grid, the nearest
   point within the grid is still returned.

   Arguments:
      lats (np.array): 2D array of latitudes
      lons (np.array): 2D array of longitude
      lat (float): Loopup latitude
      lon (float): Loopup longitude

   Returns:
      I (int): First index into lats/lons arrays
      J (int): Second index into lats/lons arrays
   """
   dist = distance(lat, lon, lats, lons)
   indices = np.unravel_index(dist.argmin(), dist.shape)
   X = lats.shape[0]
   Y = lats.shape[1]
   I = indices[0]
   J = indices[1]
   if(indices[0] == 0 or indices[0] >= X-1 or indices[1] == 0 or indices[1] >= Y-1):
      warning("Lat/lon %g,%g outside grid" % (lat, lon))
   return I, J


def distance(lat1, lon1, lat2, lon2):
   """
   Computes the great circle distance between two points using the
   haversine formula. Values can be vectors.
   """
   # Convert from degrees to radians
   pi = 3.14159265
   lon1 = lon1 * 2 * pi / 360
   lat1 = lat1 * 2 * pi / 360
   lon2 = lon2 * 2 * pi / 360
   lat2 = lat2 * 2 * pi / 360
   dlon = lon2 - lon1
   dlat = lat2 - lat1
   a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
   c = 2 * np.arcsin(np.sqrt(a))
   distance = 6.367e6 * c
   return distance
