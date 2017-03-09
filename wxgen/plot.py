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
      self.which_vars = None

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
            # Keep tany set ticks and labels
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
      sim = sims[0]
      for m in range(sim.num):
         traj = sim.get(m)
         values = sim.extract(traj)
         for i in range(values.shape[1]):
            mpl.subplot(values.shape[1], 1,i+1)
            mpl.plot(values[:,i])
            mpl.ylabel(sim.variables[i].name)

      traj = truth.get(0)
      values = truth.extract(traj)
      for i in range(values.shape[1]):
         mpl.plot(values[:,i], lw=5, color="red")
         mpl.subplot(values.shape[1], 1,i+1)
      self._finish_plot()


class Variance(Plot):
   def plot(self, sims, truth):
      Ivar = 0
      scales = [1,3,7,11,31,61, 181, 365]
      if truth is not None:
         traj = truth.get(0)
         values = truth.extract(traj)[:,Ivar]
         truth_var = self.compute_truth_variance(values, scales)
         mpl.plot(scales, truth_var, 'o-')
         mpl.title(truth.variables[Ivar].name)
         mpl.ylabel("Variance (%s)" % truth.variables[Ivar].units)

      if sims is not None:
         sim = sims[0]
         sim_values = np.zeros([sim.length, sim.num])
         for m in range(sim.num):
            traj = sim.get(m)
            q = sim.extract(traj)
            sim_values[:,m] = q[:,Ivar]
         sim_var = self.compute_sim_variance(sim_values, scales)
         mpl.plot(scales, sim_var, 'o-')
      mpl.xlabel("Time scale (days)")
      mpl.grid()
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
      for i in range(0, N):
         truth[:,i] = truth[:,i] - clim

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
      for i in range(0, N):
         array[:,i] = array[:,i] - clim

      # Compute variance
      variance = np.zeros(len(scales))
      for i in range(0, len(scales)):
         s = scales[i]
         c = [1.0/s]* s
         sim_c = np.zeros([array.shape[0], N], float)
         for e in range(0, N):
            sim_c[:,e] = astropy.convolution.convolve(array[:,e], 1.0/s*np.ones(s))
         variance[i] = np.nanvar(sim_c)
      return variance
