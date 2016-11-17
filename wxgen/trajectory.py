import numpy as np
import wxgen.database


class Trajectory(object):
   """
   Represents a trajectory of states. This is represented as a sequence of indicies into the
   database: (Ensemble member index, day index).

   The trajectory can be extracted either as an aggregated sequence or the full gridded sequence.
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

   def extract(self):
      """ Returns the sequence as a 2D numpy array (T, V) """
      T = self.indices.shape[0]
      V = len(self.database.variables)
      trajectory = np.zeros([T, V], float)
      for i in range(0, self.indices.shape[0]):
         trajectory[i,:] = self.database._data_agg[self.indices[i,1],:,self.indices[i,0]]
      return trajectory

   def extract_gridded(self):
      """ Returns the sequence as a 4D numpy array (T, V, X, Y) """
      pass
