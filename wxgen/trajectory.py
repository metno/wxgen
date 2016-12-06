import numpy as np
import wxgen.database


class Trajectory(object):
   """
   Represents a trajectory of states. This is represented as a sequence of indicies into the
   database: (Ensemble member index, day index).

   The trajectory can be extracted either as an aggregated sequence or the full gridded sequence.

   Attributes:
   indices        Indices for database
   length         Length of trajectory
   variables      Variables
   """
   def __init__(self, indices, database):
      """
      Arguments:
      indices     A two column numpy array of ints. The first column is the trajectory index and the
                  second is the day index.
      database    Of type wxgen.database.Database
      """
      assert(len(indices.shape) == 2)
      assert(indices.shape[1] == 2)
      self.indices = indices
      self.database = database

   @property
   def length(self):
      return self.indices.shape[0]

   @property
   def variables(self):
      return self.database.variables

   @property
   def X(self):
      return self.database.X

   @property
   def Y(self):
      return self.database.Y

   def extract(self):
      """ Returns the sequence as a 2D numpy array (T, V) """
      T = self.indices.shape[0]
      V = len(self.database.variables)
      trajectory = np.nan*np.zeros([T, V], float)
      for i in range(0, self.indices.shape[0]):
         if not np.isnan(self.indices[i,1]):
            trajectory[i,:] = self.database._data_agg[self.indices[i,1],:,self.indices[i,0]]
      return trajectory

   def extract_grid(self):
      """ Returns the sequence as a 4D numpy array (T, X, Y, V) """
      T = self.indices.shape[0]
      V = len(self.database.variables)
      X = self.database.X
      Y = self.database.Y
      trajectory = np.zeros([T, X, Y, V], float)
      # Loop over member, lead-time indices
      for i in range(0, self.indices.shape[0]):
         m = self.indices[i,0]
         t = self.indices[i,1]
         trajectory[i,:,:,:] = self.database._data[t, :, :, :, m]
      return trajectory

   def __str__(self):
      str = ""
      for i in range(0, self.indices.shape[0]):
         str = "%s[%d,%d]," % (str, self.indices[i,0], self.indices[i,1])
      return str
