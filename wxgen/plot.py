import inspect
import netCDF4
import matplotlib.pylab as mpl
import numpy as np
import sys
import wxgen.util
import matplotlib.dates
import scipy.ndimage
import astropy.convolution
import datetime
import copy


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

class Plot(object):
   def __init__(self):
      self.filename = None
      self.dpi = 300
      self.fig_size = [10,5]
      self.xlim = None
      self.ylim = None
      self.xlog = False
      self.ylog = False
      self.xticks = None
      self.yticks = None
      self.xticklabels = None
      self.yticklabels = None
      self._sets_xticks = None
      self._sets_yticks = None
      self.vars = None
      self.thresholds = None

   def plot(self, sims, truth):
      """
      Arguments:
         sims (list): List of simulation databases
         truth (wxgen.database): Truth database
      """
      raise NotImplementedError()

   def _finish_plot(self):
      for ax in mpl.gcf().get_axes():
         if self.xlim is not None:
            ax.set_xlim(self.xlim)
         if self.ylim is not None:
            ax.set_ylim(self.ylim)
         if self.xlog:
            # Keep any set ticks and labels
            if self._sets_xticks:
               xticks = ax.get_xticks()
               xticklabels = ax.get_xticklabels()
            ax.set_xscale("log")
            if self._sets_xticks:
               ax.set_xticks(xticks)
               ax.set_xticklabels(xticklabels)
         if self.ylog:
            if self._sets_yticks:
               yticks = ax.get_yticks()
               yticklabels = ax.get_yticklabels()
            ax.set_yscale("log")
            if self._sets_yticks:
               ax.set_yticks(yticks)
               ax.set_yticklabels(yticklabels)
      if self.filename is None:
         mpl.show()
      else:
         mpl.gcf().set_size_inches(self.fig_size[0],self.fig_size[1])
         mpl.savefig(self.filename, bbox_inches='tight', dpi=self.dpi)


class Timeseries(Plot):
   def plot(self, sims, truth):
      if self.vars is None:
         if truth is not None:
            Ivars = range(len(truth.variables))
         else:
            Ivars = range(len(sims[0].variables))
      else:
         Ivars = self.vars

      X = (truth is not None)
      if sims is not None:
         X += len(sims)
      Y = len(Ivars)
      if sims is not None:
         for s in range(len(sims)):
            sim = sims[s]
            for m in range(sim.num):
               traj = sim.get(m)
               values = sim.extract(traj)
               for i in range(len(Ivars)):
                  index = s*Y+i+1
                  mpl.subplot(X,Y,index)
                  Ivar = Ivars[i]
                  mpl.plot(values[:,Ivar], 'o-')
                  mpl.ylabel(sim.variables[Ivar].name)
                  mpl.title(index)

      if truth is not None:
         traj = truth.get(0)
         values = truth.extract(traj)
         for i in range(len(Ivars)):
            index = (X-1)*Y+i+1
            mpl.subplot(X,Y,index)
            Ivar = Ivars[i]
            mpl.plot(values[:,Ivar], lw=5, color="red")
            mpl.ylabel(truth.variables[Ivar].name)

      self._finish_plot()


class Variance(Plot):
   def __init__(self):
      Plot. __init__(self)
      self._sets_xticks = True
      self._normalize = False

   def plot(self, sims, truth):
      if self.thresholds is None:
         scales = [1,3,7,11,31,61, 181, 365]
      else:
         scales = self.thresholds

      if self.vars is None:
         if truth is not None:
            Ivars = range(len(truth.variables))
         else:
            Ivars = range(len(sims[0].variables))
      else:
         Ivars = self.vars
      for i in range(len(Ivars)):
         Ivar = Ivars[i]
         mpl.subplot(1, len(Ivars), i+1)
         if truth is not None:
            traj = truth.get(0)
            values = truth.extract(traj)[:,Ivar]
            # I = np.where(np.isnan(values) == 0)[0]
            # values = values[I]
            truth_var = self.compute_truth_variance(values, scales)
            mpl.plot(scales, truth_var, 'o-', label="Truth")
            mpl.title(truth.variables[Ivar].name)
            if self._normalize:
               mpl.ylabel("Normalized variance")
            else:
               mpl.ylabel("Variance ($%s^s$)" % truth.variables[Ivar].units)

         if sims is not None:
            for sim in sims:
               sim_values = np.zeros([sim.length, sim.num])
               for m in range(sim.num):
                  traj = sim.get(m)
                  q = sim.extract(traj)
                  sim_values[:,m] = q[:,Ivar]
               sim_var = self.compute_sim_variance(sim_values, scales)
               mpl.plot(scales, sim_var, 'o-', label=sim.name)
         ticks = np.array([1,7,30,365])
         labels = ["day", "week", "month", "year"]
         I = np.where(ticks < mpl.xlim()[1])[0]
         # Include the first one above the limits
         if len(I) < len(ticks):
           I = np.append(I, I[-1]+1)

         mpl.gca().set_xticks(ticks[I])
         mpl.gca().set_xticklabels(labels)
         mpl.xlabel("Time scale (days)")
         mpl.grid()
      mpl.legend()
      self._finish_plot()

   def compute_truth_variance(self, array, scales):
      """
         array: 1D array
      """
      # Create 1-year long segments
      N = int(np.ceil(len(array)/365))
      truth = np.zeros([365, N])
      for i in range(0, N):
         I = range(i*365, (i+1)*365)
         truth[:,i] = array[I]

      # Remove climatology so we can look at annomalies. Use separate obs and fcst climatology
      # otherwise the fcst variance is higher because obs gets the advantage of using its own
      # climatology.
      clim = np.nanmean(truth, axis=1)
      std = 1
      if self._normalize:
         std = np.nanstd(truth, axis=1)

      for i in range(0, N):
         truth[:,i] = (truth[:,i] - clim)/std

      # Compute variance
      variance = np.zeros(len(scales))
      for i in range(0, len(scales)):
         s = scales[i]
         c = [1.0/s]* s
         truth_c = np.zeros([truth.shape[0], N], float)
         for e in range(0, N):
            truth_c[:,e] = astropy.convolution.convolve(truth[:,e], 1.0/s*np.ones(s))
         variance[i] = np.nanvar(truth_c)
      return variance

   def compute_sim_variance(self, array, scales):
      """
      Arguments:
         array (np.array): 2D array (time, member)
         scales (list): List of time lengths

      Returns:
         list: Variance for different time lengths
      """
      # Create 1-year long segments
      N = array.shape[1]

      # Remove climatology so we can look at annomalies. Use separate obs and fcst climatology
      # otherwise the fcst variance is higher because obs gets the advantage of using its own
      # climatology.
      clim = np.nanmean(array, axis=1)
      std = 1
      if self._normalize:
         std = np.nanstd(array, axis=1)
      values = copy.deepcopy(array)
      for i in range(0, N):
         values[:,i] = (values[:,i] - clim)/std

      # Compute variance
      variance = np.nan*np.zeros(len(scales))
      for i in range(0, len(scales)):
         s = scales[i]
         if array.shape[0] > s:
            c = [1.0/s]* s
            sim_c = np.zeros([values.shape[0], N], float)
            for e in range(0, N):
               sim_c[:,e] = astropy.convolution.convolve(values[:,e], 1.0/s*np.ones(s))
            variance[i] = np.nanvar(sim_c)
      return variance
