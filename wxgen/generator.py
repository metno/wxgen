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
      Returns a list of N trajectories, where each trahectry has a length of T.

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
      X = self._database.X
      Y = self._database.Y

      for n in range(0, N):
         trajectory_indices = np.zeros([T, 2], int)

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
            month_of_year = day_of_year / 30

            segment_curr = self._database.get_random(state_curr, self._metric, month_of_year)
            indices_curr = segment_curr.indices

            end = min(start + Tsegment-1, T)  # Ending index
            Iout = range(start, end)  # Index into trajectory
            Iin = range(0, end - start)  # Index into segment
            trajectory_indices[Iout, :] = indices_curr[Iin, :]
            if self._debug:
               print "Current state: ", state_curr
               print "Chosen segment: ", segment_curr
               print "Trajectory indices: ", Iout
               print "Segment indices: ", Iin
            state_curr = segment_curr.extract()[-1,:]
            start = start + Tsegment-1
            day_of_year = day_of_year + Tsegment-1
            if day_of_year > 365:
               day_of_year = day_of_year - 365

         trajectory = wxgen.trajectory.Trajectory(trajectory_indices, self._database)
         if self._debug:
            print "Trajectory: ", trajectory
         trajectories.append(trajectory)

      return trajectories
