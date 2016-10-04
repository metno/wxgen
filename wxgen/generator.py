import numpy as np
import wxgen.metric
import wxgen.util
class Generator(object):
   def __init__(self, database):
      self._database = database

   #@profile
   def get(self, N, T, initial_state=None):
      # Initialize
      trajectories = list()
      V = self._database.num_vars()
      for n in range(0, N):
         trajectory = np.zeros([T, V], float)

         T0 = self._database.days()
         assert(T % T0 == 0)

         # Assemble a trajectory by concatenating appropriate trajectories
         counter = 0
         if initial_state is None:
            state_curr = self._database.get(0)[0,:]
         else:
            state_curr = initial_state

         while counter < T/T0:
            segment_curr = self._database.get_random(state_curr, wxgen.metric.Rmsd)

            indices = range(counter*T0, (counter+1) * T0)
            for var in self._database.vars():
               trajectory[indices, :] = segment_curr
            counter = counter + 1
            state_curr = segment_curr[-1, :]

         trajectories.append(trajectory)
      return trajectories
