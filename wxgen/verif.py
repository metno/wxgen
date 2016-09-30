import numpy as np

# Compute verification statistics on trajectories
class Verif:
   pass

class Change(Verif):
   def compute(self, trajectories):
      var = "T"
      N = len(trajectories)
      T = trajectories[0][var].shape[0]
      changes = np.zeros(T-1, float)
      for t in range(0, T-1):
         ar = np.array([abs(trajectories[i][var][t] - trajectories[i][var][t+1]) for i in range(0, N)])
         changes[t] = np.mean(ar)
      return changes

