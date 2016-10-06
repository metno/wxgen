import inspect
import matplotlib.pylab as mpl
import numpy as np
import sys
import wxgen.util


def get_all():
   """ Returns a list of all output classes """
   temp = inspect.getmembers(sys.modules[__name__], inspect.isclass)
   return temp


def get(name):
   """ Returns an output object of a class with the given name """
   outputs = get_all()
   m = None
   for mm in outputs:
      if(name == mm[0].lower()):
         m = mm[1]
   if m is None:
      wxgen.util.error("Cannot find output called '%s'" % name)
   return m


class Output(object):
   """
   A class for outputing trajectory information
   """
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
   """
   Draws all trajectories as lines. One variable per subplot.
   """
   def plot(self, trajectories):
      T = trajectories[0].shape[0]
      V = trajectories[0].shape[1]
      vars = self._db.vars()
      Tsegment = self._db.days()
      for v in range(0, V):
         mpl.subplot(V,1,v+1)
         x = np.linspace(0, T-1, T)
         for tr in trajectories:
            # Plot the trajectory
            mpl.plot(x, tr[:,v], 'k.-', lw=0.5)

            # Plot the starting state of each segment
            I = range(0, T, Tsegment-1)
            mpl.plot(x[I], tr[I,v], 'ko', mfc='w')

         mpl.xlabel("Time (days)")
         mpl.ylabel(vars[v])
         mpl.grid()
         mpl.xlim([0, T])
      mpl.xlabel("Time (days)")
      self._finish_plot()


class Text(Output):
   """
   Writes the trajectories to a text file. One variable in each column and each day on a separate
   line. Trajectories are separated by a blank line.
   """
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


class Verification(Output):
   """
   Plots verification data for the trajectories.
   """
   _pool = True
   def plot(self, trajectories):
      N = len(trajectories)
      T = trajectories[0].shape[0]
      V = trajectories[0].shape[1]
      Tsegment = self._db.days()
      vars = self._db.vars()
      if self._pool:
         size = Tsegment-1
         x = np.linspace(0, Tsegment-2, Tsegment-1)
      else:
         size = T
         x = np.linspace(0, T-1, T)
      changes = np.zeros(size, float)
      counter = np.zeros(size, int)
      for v in range(0, V):
         mpl.subplot(V, 1, v+1)
         mpl.title(vars[v])
         for t in range(0, T-1):
            ar = np.array([abs(trajectories[i][t,v] - trajectories[i][t+1,v]) for i in range(0, N)])
            if self._pool:
               I = t % (Tsegment-1)  # use this to pool similar leadtimes
            else:
               I = t
            changes[I] = changes[I] + np.mean(ar)
            counter[I] = counter[I] + 1
         mpl.plot(x, changes / counter, 'k-')
         mpl.xlabel("Day")
         mpl.ylabel("Average absolute change to next day")
         mpl.grid()
         mpl.xlim([0, size])
         ylim = mpl.ylim()
         mpl.ylim([0, ylim[1]])
      self._finish_plot()

