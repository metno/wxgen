import numpy as np
import wxgen.metric
import wxgen.util
class Generator(object):
   def __init__(self, database, metric=wxgen.metric.Rmsd()):
      self._database = database
      self._metric = metric

   #@profile
   def get(self, N, T, initial_state=None):
      # Initialize
      trajectories = list()
      V = self._database.num_vars()
      Tsegment = self._database.days()

      for n in range(0, N):
         trajectory = np.zeros([T, V], float)

         counter = 0
         if initial_state is None:
            state_curr = self._database.get(0)[0,:]
         else:
            state_curr = initial_state

         # Assemble a trajectory by concatenating appropriate segments. Start with the
         # initial state, then find a segment that has a similar starting state. Discard the
         # initial state of this segment, and insert the remaining part of the segment.
         # Repeat. This means that if the segment is 10 days long, we are only using 9 days
         # of the segment.
         start = 1  # Starting index into output trajectory where we are inserting a segment
         trajectory[0, :] = state_curr
         while start < T:
            segment_curr = self._database.get_random(state_curr, self._metric)

            end = min(start + Tsegment-1, T)  # Ending index
            Iout = range(start, end)  # Index into trajectory
            Iin = range(1, end - start+1)  # Index into segment
            trajectory[Iout, :] = segment_curr[Iin, :]
            state_curr = segment_curr[-1, :]
            start = start + Tsegment-1

         trajectories.append(trajectory)
      return trajectories
