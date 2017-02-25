import numpy as np
import wxgen.database


class Trajectory(object):
   """
   Represents a trajectory of states. This is represented as a sequence of indicies into some
   database: (Ensemble member index, day index).

   The trajectory can be extracted either as an aggregated sequence or the full gridded sequence.

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

   def extract(self, database):
      """ Returns the sequence as a 2D numpy array (T, V) """
      T = self.indices.shape[0]
      V = len(database.variables)
      trajectory = np.nan*np.zeros([T, V], float)
      for i in range(0, self.indices.shape[0]):
         if self.indices[i,1] >= 0:
            trajectory[i,:] = database._data_agg[self.indices[i,1],:,self.indices[i,0]]
      return trajectory

   def extract_grid(self, database):
      """ Returns the sequence as a 4D numpy array (T, X, Y, V) """
      T = self.indices.shape[0]
      V = len(database.variables)
      X = database.X
      Y = database.Y
      trajectory = np.zeros([T, X, Y, V], float)
      # Loop over member, lead-time indices
      for i in range(0, self.indices.shape[0]):
         m = self.indices[i,0]
         t = self.indices[i,1]
         assert(not np.isnan(m))
         assert(not np.isnan(t))
         trajectory[i,:,:,:] = database._data[t, :, :, :, m]
      return trajectory

   def __str__(self):
      str = ""
      for i in range(0, self.indices.shape[0]):
         str = "%s[%d,%d]," % (str, self.indices[i,0], self.indices[i,1])
      return str
