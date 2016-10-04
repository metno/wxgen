import numpy as np
import wxgen.metric
import wxgen.util
class Generator(object):
   def __init__(self, database, metric=wxgen.metric.Rmsd()):
      self._database = database
      self._metric = metric

   #@profile
   def get(self, N, T, initial_state=None):
      debug = False
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

         # Assemble a trajectory by concatenating appropriate segments. Start by finding a
         # segment that has a starting state that is similar to the requested initial state.
         # When repeating, overwrite the end state of the previous segment. This means that
         # if the segment is 10 days long, we are only using 9 days of the segment.
         start = 0  # Starting index into output trajectory where we are inserting a segment
         while start < T:
            segment_curr = self._database.get_random(state_curr, self._metric)

            end = min(start + Tsegment-1, T)  # Ending index
            Iout = range(start, end)  # Index into trajectory
            Iin = range(0, end - start)  # Index into segment
            trajectory[Iout, :] = segment_curr[Iin, :]
            if debug:
               print "Current state: ", state_curr
               print "Chosen segment: ", segment_curr
               print "Trajectory indices: ", Iout
               print "Segment indices: ", Iin
            state_curr = segment_curr[-1, :]
            start = start + Tsegment-1

         if debug:
            print "Trajectory: ", trajectory
         trajectories.append(trajectory)
      return trajectories
