import matplotlib.pylab as mpl
import wxgen.util


class Output(object):
   def __init__(self, filename=None):
      self._filename = filename
      self._dpi = 300

   def _finish_plot(self):
      if self._filename is None:
         mpl.show()
      else:
         mpl.savefig(self._filename, dpi=self._dpi)


class Timeseries(Output):
   def plot(self, trajectories):
      T = trajectories[0].shape[0]
      v = 0
      for tr in trajectories:
         mpl.plot(tr[:,v], 'k.-', lw=0.5)
      mpl.xlabel("Time (days)")
      mpl.ylabel("Value")
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
