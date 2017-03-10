import numpy as np
import calendar
import datetime
import sys


DEBUG = False


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


def debug(message):
   if DEBUG:
      print "\033[1;32mDebug: " + message + "\033[0m"


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
               values = values + list(np.arange(start, end+ stepSign*0.0001, step))
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


