import numpy as np
import sys

# Randomly selects an index in an array based on weights
def random_weighted(weights):
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


# allowable formats:
# num
# num1,num2,num3
# start:end
# start:step:end
def parse_numbers(numbers, isDate=False):
   # Check if valid string
   if(any(char not in set('-01234567890.:,') for char in numbers)):
      error("Could not translate '" + numbers + "' into numbers")

   values = list()
   commaLists = numbers.split(',')
   for commaList in commaLists:
      colonList = commaList.split(':')
      if(len(colonList) == 1):
         values.append(float(colonList[0]))
      elif(len(colonList) <= 3):
         start = float(colonList[0])
         step = 1
         if(len(colonList) == 3):
            step = float(colonList[1])
         stepSign = step / abs(step)
         # arange does not include the end point:
         end = float(colonList[-1]) + stepSign * 0.0001
         if(isDate):
            date = min(start, end)
            curr = list()
            while date <= max(start, end):
               curr.append(date)
               date = getDate(date, step)
            values = values + list(curr)
         else:
            values = values + list(np.arange(start, end, step))
      else:
         error("Could not translate '" + numbers + "' into numbers")
      if(isDate):
         for i in range(0, len(values)):
            values[i] = int(values[i])
   return values


def resize(vec, size):
   if vec.shape[0] == size[0] and len(vec.shape) == 1:
      vec_resized = np.reshape(np.repeat(vec, size[1]), size)
   elif vec.shape[0] == 1 and len(vec.shape) == 1:
      vec_resized = vec*np.ones(size)
   else:
      # Check that the output dims are multiples of input dims
      assert(size[0] % vec.shape[0] == 0)
      assert(size[1] % vec.shape[1] == 0)
      vec_resized = np.tile(vec, (size[0] / vec.shape[0], size[1] / vec.shape[1]))
   return vec_resized
