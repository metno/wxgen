import numpy as np
import wxgen.metric
import wxgen.util
import wxgen.climate_model


class LargeScale(object):
   """
   This class generates long trajectories from segments in a database
   """

   def __init__(self, database, metric=wxgen.metric.Rmsd(), model=wxgen.climate_model.Bin(10)):
      self._database = database
      self._metric = metric
      self._model = model
      self.prejoin = None

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
            tr = self.get_random(np.zeros(V), wxgen.metric.Exp(np.zeros(V)))
            state_curr = self._database.extract(tr)[0,:]
         else:
            state_curr = initial_state

         # Assemble a trajectory by concatenating appropriate segments. Start by finding a segment
         # that has a starting state that is similar to the requested initial state. When
         # repeating, overwrite the end state of the previous segment. This means that if the
         # segment is 10 days long, we are only using 9 days of the segment.
         start = 0  # Starting index into output trajectory where we are inserting a segment
         time = wxgen.util.date_to_unixtime(20170101)
         join = 0
         while start < T:
            # print "Join: %d" % join
            # TODO
            climate_state = self._model.get([time])[0]

            search_times = None
            if join > 0:
               print 1
               end_times = self._database.inittimes[segment_curr.indices[-1,0]]
               search_times = [end_times - 5*86400, end_times + 5*86400]
            segment_curr = self.get_random(state_curr, self._metric, climate_state, search_times)
            indices_curr = segment_curr.indices

            end = min(start + Tsegment-1, T)  # Ending index
            Iout = range(start, end)  # Index into trajectory
            Iin = range(0, end - start)  # Index into segment
            trajectory_indices[Iout, :] = indices_curr[Iin, :]
            wxgen.util.debug("Current state: %s" % state_curr)
            wxgen.util.debug("Chosen segment: %s" % segment_curr)
            wxgen.util.debug("Trajectory indices: %s" % Iout)
            wxgen.util.debug("Segment indices: %s" % Iin)
            state_curr = self._database.extract(segment_curr)[-1,:]
            start = start + Tsegment-1
            time = time + (Tsegment-1)*86400
            if self.prejoin > 0:
               join = (join + 1) % self.prejoin

         trajectory = wxgen.trajectory.Trajectory(trajectory_indices)
         wxgen.util.debug("Trajectory: %s" % trajectory)
         trajectories.append(trajectory)

      return trajectories

   def get_random(self, target_state, metric, climate_state=None, time_range=None):
      """
      Returns a pseudo-random segment from the database chosen based on weights computed by a metric

      Arguments:
         target_state (np.array): Try to matchin this state when finding the trajectory. One value
            for each variable in the database.
         metric (wxgen.metric): Metric to use when finding matches
         climate_state (np.array): External state representing what state the climate is in
         time_range (list): Start and end unixtimes for the search

      Returns:
         wxgen.trajectory: Random trajectory
      """
      assert(np.sum(np.isnan(target_state)) == 0)
      weights = metric.compute(target_state, self._database._data_agg[0,:,:])
      if time_range is None:
         Ivalid = np.where(np.isnan(weights) == 0)[0]
      else:
         Ivalid = np.where((np.isnan(weights) == 0) & (self._database.inittimes > time_range[0]) &
               (self._database.inittimes < time_range[1]))[0]
         # print Ivalid
         # print "Restricting time range"
      if climate_state is not None and time_range is None:
         Iclimate_state = np.where(self._database.climate_states[Ivalid] == climate_state)[0]
         if len(Iclimate_state) == 0:
            wxgen.util.error("Cannot find a segment with climate state = %s" % str(climate_state))
         Ivalid = Ivalid[Iclimate_state]

      weights_v = weights[Ivalid]

      # Flip the metric if it is negative oriented
      if metric._orientation == -1:
         I0 = np.where(weights_v < 1e-3)[0]
         I1 = np.where(weights_v >= 1e-3)[0]
         # Ensure we do not get too high weights
         weights_v[I1] = 1.0/weights_v[I1]
         weights_v[I0] = 1e3

      # max_weight = np.max(weights_v)
      # weights_v[np.where(weights_v < max_weight / 4)[0]] = 0
      I_v = wxgen.util.random_weighted(weights_v)
      I = Ivalid[I_v]

      # Do a weighted random choice of the weights
      wxgen.util.debug("I: %s" % I)
      wxgen.util.debug("Data: %s" % self._database._data_agg[0,:,I])
      wxgen.util.debug("Weight: %s" % weights_v[I_v])
      wxgen.util.debug("Max weight: %s" % np.max(weights_v))
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
