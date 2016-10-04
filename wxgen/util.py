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
