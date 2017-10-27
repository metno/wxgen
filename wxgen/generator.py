import numpy as np
import wxgen.metric
import wxgen.util
import wxgen.climate_model


class LargeScale(object):
   """
   This class generates long trajectories from segments in a database
   """

   def __init__(self, database, metric=wxgen.metric.Rmsd(), model=None):
      self._database = database
      self._metric = metric
      self.prejoin = None
      self.policy = "random"
      self.stagger = False

   def get(self, N, T, initial_state=None):
      """
      Returns a list of N trajectories, where each trajectory has a length of T.

      If initial_state is provided then the trajectory will start with a state similar to this. If
      no initial state is provided, start with a random segment from the database.

      Arguments:
         N (int): Number of trajectories
         T (int): Number of timesteps in each trajectory
         initial_state (np.array): An array of the initial state (must be of length V)
      """
      trajectories = list()
      V = len(self._database.variables)
      if self._database.length < 2:
         wxgen.util.error("Cannot create simulation with a database with segments shorter than 2 days")
      X = self._database.X
      Y = self._database.Y

      for n in range(0, N):
         wxgen.util.debug("Generating trajectory %d/%d" % (n+1, N), color="red")
         trajectory_indices = -1+np.zeros([T, 2], int)

         time = wxgen.util.date_to_unixtime(20170101)
         climate_state = self._database.model.get([time])[0]
         if initial_state is None:
            wxgen.util.debug("Finding random starting state", color="yellow")
            I = np.random.randint(self._database.num)
            num_vars = self._database._data_matching.shape[1]
            tr = self.get_random(np.zeros(num_vars), wxgen.metric.Exp(np.zeros(num_vars)), climate_state)
            if 1:
               state_curr = self._database.extract_matching(tr)[0, :]
            else:
               state_curr = self._database.extract(tr)[0, :]
         else:
            state_curr = initial_state

         """
         Assemble a trajectory by concatenating appropriate segments. Start by finding a segment
         that has a starting state that is similar to the requested initial state. When repeating,
         overwrite the end state of the previous segment. This means that if the segment is 10 days
         long, we are only using 9 days of the segment.

         If the last segment fits exactly, we do not need to search for a new state to fill in the
         last day with new data. This is the case when start == T-1.
         """
         start = 0  # Starting index into output trajectory where we are inserting a segment
         join = 0
         while start < T-1:
            climate_state = self._database.model.get([time])[0]

            """
            Prejoin multiple segments that are nearby in time. This is done by passing
            'search_times' to get_random.
            """
            search_times = None
            if join > 0:
               end_times = self._database.inittimes[segment_curr.indices[-1, 0]] + segment_curr.indices[-1, 1]*86400
               search_times = [end_times - 5*86400, end_times + 5*86400]
            wxgen.util.debug("Found random segment", color="yellow")
            wxgen.util.debug("Target state: %s" % ' '.join(["%0.2f" % x for x in state_curr]))
            segment_curr = self.get_random(state_curr, self._metric, climate_state, search_times)

            """
            Stagger the trajectories so that they don't all potentially have jumps at the same
            leadtimes. This is done by truncating the first segment to a random length. Note that
            the upper end of randint is exclusive, hence the "+ 1".
            """
            if self.stagger and start == 0:
               I = np.random.randint(1, segment_curr.indices.shape[0] + 1)
               segment_curr.indices = segment_curr.indices[0:I, :]

            indices_curr = segment_curr.indices

            """
            Account for the fact that the desired trajectory length is not a whole multiple of the
            segment length: Only take the first part of the segment if needed.
            """
            Tsegment = len(indices_curr)
            end = min(start + Tsegment, T)  # Ending index
            Iout = range(start, end)  # Index into trajectory
            Iin = range(0, end - start)  # Index into segment
            trajectory_indices[Iout, :] = indices_curr[Iin, :]

            # wxgen.util.debug("Current state: %s" % state_curr)
            # wxgen.util.debug("Chosen segment: %s" % segment_curr)
            # wxgen.util.debug("Trajectory indices: %s" % Iout)
            # wxgen.util.debug("Segment indices: %s" % Iin)
            state_curr = self._database.extract_matching(segment_curr)[-1, :]
            start = start + Tsegment-1
            time = time + (Tsegment-1)*86400
            if self.prejoin is not None and self.prejoin > 0:
               join = (join + 1) % self.prejoin

         if len(np.where(trajectory_indices == -1)[0]) > 0:
            wxgen.util.error("Internal error. The trajectory was not properly filled")
         trajectory = wxgen.trajectory.Trajectory(trajectory_indices)
         wxgen.util.debug("Trajectory: %s" % trajectory)
         trajectories.append(trajectory)

      return trajectories

   def get_random(self, target_state, metric, climate_state=None, time_range=None):
      """
      Returns a pseudo-random segment from the database chosen based on weights computed by a metric

      Arguments:
         target_state (np.array): Try to match this state when finding the trajectory. One value
            for each variable in the database.
         metric (wxgen.metric): Metric to use when finding matches
         climate_state (np.array): External state representing what state the climate is in
         time_range (list): Start and end unixtimes for the search. If None, then do not restrict.

      Returns:
         wxgen.trajectory: Random trajectory
      """
      assert(np.sum(np.isnan(target_state)) == 0)

      weights = metric.compute(target_state, self._database._data_matching[0, :, :])
      use_climate_state = climate_state is not None

      # Find valid segments
      do_prejoin = False
      if time_range is None:
         Itime = np.where(np.isnan(weights) == 0)[0]
      else:
         do_prejoin = True
         Itime = np.where((np.isnan(weights) == 0) & (self._database.inittimes > time_range[0]) & (self._database.inittimes < time_range[1]))[0]
         if len(Itime) == 0:
            date_range = [wxgen.util.unixtime_to_date(t) for t in time_range]
            wxgen.util.warning("Skipping this prejoin: No valid segment that start in date range [%d, %d]" %
                  (date_range[0], date_range[1]))
            Itime = np.where(np.isnan(weights) == 0)[0]
            # Without any prejoin segments, revert to the original plan of just finding a random segment
            do_prejoin = False

      # Segment the database based on climate state
      if climate_state is not None and not do_prejoin:
         Iclimate_state = np.where(self._database.climate_states[Itime] == climate_state)[0]
         if len(Iclimate_state) == 0:
            wxgen.util.error("Cannot find a segment with climate state = %s" % str(climate_state))
         Itime = Itime[Iclimate_state]

      weights_v = weights[Itime]

      # Flip the metric if it is negative oriented
      if metric._orientation == -1:
         I0 = np.where(weights_v < 1e-3)[0]
         I1 = np.where(weights_v >= 1e-3)[0]
         # Ensure we do not get too high weights
         weights_v[I1] = 1.0/weights_v[I1]
         weights_v[I0] = 1e3

      I_v = wxgen.util.random_weighted(weights_v, self.policy)
      I = Itime[I_v]

      # Do a weighted random choice of the weights
      wxgen.util.debug("Found state:  %s" % ' '.join(["%0.2f" % x for x in self._database._data_matching[0, :, I]]))
      wxgen.util.debug("Date: %s (%i)" % (wxgen.util.unixtime_to_date(self._database.inittimes[I]), I))
      wxgen.util.debug("Climate: %s" % (climate_state))
      wxgen.util.debug("Weight (max weight): %s (%s)" % (weights_v[I_v], np.max(weights_v)))
      return self._database.get(I)


class SmallScale(object):
   """
   Generates high-resolution gridded fields from a large-scale trajectory
   """
   def __init__(self, downscalers):
      self.downscalers = downscalers

   def extract(self, trajectory, database):
      """
      Creates a high resolution gridded array

      Returns:
         data (np.array): A 4D array (T, X, Y, V)
      """
      values = self.downscale(trajectory, database)

      for downscaler in self.downscalers:
         values = downscaler.generate(values)

      return values

   def write(self, output):
      raise NotImplementedError()

   def downscale(self, trajectory, database, grid):
      """
      Downscale grid to higher resolution
      """
      large_scale = database.extract_grid(trajectory)
      raise NotImplementedError()
