import numpy as np
import metric
import util
class Generator:
   def __init__(self, database):
      self._database = database

   #@profile
   def get(self, N, T, initial_state=None):
      # Initialize
      trajectories = list()
      for n in range(0, N):
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
         for var in curr:
            curr[var] = np.zeros(1, float)+curr[var][0]
         weights = np.zeros(len(db), float)
         while counter < T/T0:
            # Find best match
            state_end = dict()
            for var in self._database.vars():
               state_end[var] = curr[var][-1]

            indices = range(counter*T0, (counter+1) * T0)
            curr = self._database.get_min(state_end, metric.Rmsd)
            for var in self._database.vars():
               trajectory[var][indices] = curr[var]
            counter = counter + 1

         trajectories.append(trajectory)
      return trajectories
