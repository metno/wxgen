import numpy as np


class Trajectory(object):
   """
   Represents a trajectory of states. This is represented as a sequence of indicies into some
   database: (Ensemble member index, day index).

   Attributes:
      indices (np.array): Indices for database
      length (int): Length of trajectory
   """
   def __init__(self, indices):
      """
      Arguments:
         indices (np.array): Array with two columns ints. The first column is the trajectory index
            and the second is the day index.
      """
      assert(len(indices.shape) == 2)
      assert(indices.shape[1] == 2)
      self.indices = indices

   @property
   def length(self):
      return self.indices.shape[0]

   def __str__(self):
      str = ""
      for i in range(0, self.indices.shape[0]):
         str = "%s[%d,%d]," % (str, self.indices[i,0], self.indices[i,1])
      return str
