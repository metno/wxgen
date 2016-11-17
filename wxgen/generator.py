import numpy as np
import wxgen.metric
import wxgen.util


class Generator(object):
   """
   This class generates long trajectories from segments in a database
   """

   def __init__(self, database, metric=wxgen.metric.Rmsd()):
      self._database = database
      self._metric = metric
      self._debug = False

   def get(self, N, T, initial_state=None):
      """
      Returns a 2D numpy array of trajectories. Each value is an index into the date and member of
      the database.

      If initial_state is provided then the trajectory will start with a state similar to this. If
      no initial state is provided, start with a random segment from the database.

      N              Number of trajectories
      T              Number of timesteps in each trajectory
      initial_state  A numpy array of the initial state (must be of length V)
      """
      # Initialize
      trajectories = list()
      V = len(self._database.variables)
      Tsegment = self._database.length

      for n in range(0, N):
         trajectory = np.zeros([T, V], float)

         if initial_state is None:
            I = np.random.randint(self._database.num)
            state_curr = self._database.get(I)
         else:
            state_curr = initial_state

         # Assemble a trajectory by concatenating appropriate segments. Start by finding a segment
         # that has a starting state that is similar to the requested initial state. When
         # repeating, overwrite the end state of the previous segment. This means that if the
         # segment is 10 days long, we are only using 9 days of the segment.
         start = 0  # Starting index into output trajectory where we are inserting a segment
         day_of_year = 1
         while start < T:
            # TODO
            month_of_year = 0 #day_of_year / 30
            segment_curr = self._database.get_random(state_curr, self._metric, month_of_year).get()
            #segment_curr = self._database.get_random(state_curr, self._metric)

            end = min(start + Tsegment-1, T)  # Ending index
            Iout = range(start, end)  # Index into trajectory
            Iin = range(0, end - start)  # Index into segment
            trajectory[Iout, :] = segment_curr[Iin, :]
            if self._debug:
               print "Current state: ", state_curr
               print "Chosen segment: ", segment_curr
               print "Trajectory indices: ", Iout
               print "Segment indices: ", Iin
            state_curr = segment_curr[-1, :]
            start = start + Tsegment-1
            day_of_year = day_of_year + Tsegment-1
            if day_of_year > 365:
               day_of_year = day_of_year - 365

         if self._debug:
            print "Trajectory: ", trajectory
         trajectories.append(trajectory)

      return trajectories
