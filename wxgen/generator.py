import numpy as np
import metric
import util
class Generator:
   def __init__(self, database):
      self._database = database

   def get(self, T, initial_state=None):
      # Initialize
      trajectory = dict()
      for var in self._database.vars():
         trajectory[var] = np.zeros(T, float)

      T0 = self._database.days()
      assert(T % T0 == 0)

      db = [self._database.get(i) for i in range(0, self._database.size())]

      # Assemble a trajectory by concatenating appropriate trajectories
      counter = 0
      if initial_state is None:
         curr = db[0]
      else:
         curr = initial_state
      weights = np.zeros(len(db), float)
      while counter < T/T0:
         # Find best match
         state_end = dict()
         for var in self._database.vars():
            state_end[var] = curr[var][-1]
         for i in range(0, len(db)):
            state_start = dict()
            for var in self._database.vars():
               state_start[var] = db[i][var][0]
            weights[i] = 1.0/metric.Rmsd.compute(state_start, state_end)
            #print state_end, state_start, curr_score
         #print state_end, Ibest, min_score

         # Do a weighted random choice of the weights
         Inext = util.random_weighted(weights)

         indices = range(counter*T0, (counter+1) * T0)
         for var in self._database.vars():
            trajectory[var][indices] = db[Inext][var]
         counter = counter + 1
         curr = db[Inext]

      return trajectory
