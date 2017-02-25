import numpy as np
import wxgen.metric
import wxgen.util
import wxgen.climate_model


class Generator(object):
   """
   This class generates long trajectories from segments in a database
   """

   def __init__(self, database, metric=wxgen.metric.Rmsd(), model=wxgen.climate_model.Bin(10)):
      self._database = database
      self._metric = metric
      self._debug = False
      self._model = model

   def get(self, N, T, initial_state=None):
      """
      Returns a list of N trajectories, where each trahectry has a length of T.

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
            state_curr = self.get_random(np.zeros(V),
                  wxgen.metric.Exp(np.zeros(V))).extract(self._database)[0,:]
         else:
            state_curr = initial_state

         # Assemble a trajectory by concatenating appropriate segments. Start by finding a segment
         # that has a starting state that is similar to the requested initial state. When
         # repeating, overwrite the end state of the previous segment. This means that if the
         # segment is 10 days long, we are only using 9 days of the segment.
         start = 0  # Starting index into output trajectory where we are inserting a segment
         time = wxgen.util.date_to_unixtime(20170101)
         while start < T:
            # TODO
            climate_state = self._model.get([time])[0]

            segment_curr = self.get_random(state_curr, self._metric, climate_state)
            indices_curr = segment_curr.indices
            # print state_curr, segment_curr.extract(self._database)[0,:]

            end = min(start + Tsegment-1, T)  # Ending index
            Iout = range(start, end)  # Index into trajectory
            Iin = range(0, end - start)  # Index into segment
            trajectory_indices[Iout, :] = indices_curr[Iin, :]
            if self._debug:
               print "Current state: ", state_curr
               print "Chosen segment: ", segment_curr
               print "Trajectory indices: ", Iout
               print "Segment indices: ", Iin
            state_curr = segment_curr.extract(self._database)[-1,:]
            start = start + Tsegment-1
            time = time + (Tsegment-1)*86400

         trajectory = wxgen.trajectory.Trajectory(trajectory_indices)
         if self._debug:
            print "Trajectory: ", trajectory
         trajectories.append(trajectory)

      return trajectories

   def get_random(self, target_state, metric, climate_state=None):
      """
      Returns a random segment from the database that is weighted
      by the scores computed by metric.
      target_state   A numpy array (length V)
      metric         Of type wxgen.metric.Metric
      climate_state      External state

      Returns:
      trajectory     Of type wxgen.trajectory.Trajectory
      """
      weights = metric.compute(target_state, self._database._data_agg[0,:,:])
      Ivalid = np.where(np.isnan(weights) == 0)[0]
      if climate_state is not None:
         Iclimate_state = np.where(self._database.climate_states[Ivalid] == climate_state)[0]
         if len(Iclimate_state) == 0:
            print np.unique(self._database.climate_states[Ivalid])
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
      if self._database._debug:
         print "I: ", I
         print "Data: ", self._database._data_agg[0,:,I]
         print "Weight: ", weights_v[I_v]
         print "Max weight: ", np.max(weights_v)
      return self._database.get(I)

   def get_sequence(self, indices):
      """ Returns a gridded sequence of states

      Arguments:
      indices     A numpy array of integers

      Returns:
      data        A 4D array (T, X, Y, V)
      """
      pass

