import numpy as np
import matplotlib.pylab as mpl
import wxgen.util


class Output(object):
   def __init__(self, db, filename=None):
      self._filename = filename
      self._db = db
      self._dpi = 300

   def _finish_plot(self):
      if self._filename is None:
         mpl.show()
      else:
         mpl.savefig(self._filename, dpi=self._dpi)


class Timeseries(Output):
   def plot(self, trajectories):
      T = trajectories[0].shape[0]
      V = trajectories[0].shape[1]
      vars = self._db.vars()
      Tsegment = self._db.days()
      for v in range(0, V):
         mpl.subplot(V,1,v+1)
         x = np.linspace(0, T-1, T)
         for tr in trajectories:
            mpl.plot(x, tr[:,v], 'k.-', lw=0.5)
            I = range(0, T-1, Tsegment)
            mpl.plot(x[I], tr[I,v], 'ro')
         mpl.xlabel("Time (days)")
         mpl.ylabel(vars[v])
         mpl.grid()
         mpl.xlim([0, T])
         mpl.xticks(np.linspace(0, T, T / Tsegment+1))
      self._finish_plot()

class Database(Output):
   def plot(self, trajectories):
      data = self._db._data
      V = data.shape[1]
      T = data.shape[0]
      N = 1000#data.shape[2]
      vars = self._db.vars()
      for v in range(0, V):
         mpl.subplot(V,1,v+1)
         x = np.linspace(0, T, T)
         mpl.plot(x, data[:,v, :], 'k.-', lw=0.5)
         mpl.xlabel("Time (days)")
         mpl.ylabel(vars[v])
         mpl.grid()
         mpl.xlim([0, T])
      self._finish_plot()

class Text(Output):
   def plot(self, trajectories):
      if self._filename is None:
         wxgen.util.error("Text output requires a filename")

      fid = open(self._filename, "w")
      N = len(trajectories)
      T = trajectories[0].shape[0]
      V = trajectories[0].shape[1]
      for n in range(0, N):
         for t in range(0, T):
            for v in range(0, V):
               fid.write("%f " % trajectories[n][t,v])
            fid.write("\n")
         if n < N-1:
            fid.write("\n")
      fid.close()
