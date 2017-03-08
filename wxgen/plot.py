import inspect
import netCDF4
import matplotlib.pylab as mpl
import numpy as np
import sys
import wxgen.util
import matplotlib.dates
import scipy.ndimage
# import astropy.convolution
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
      traj = truth.get(0)
      values = truth.extract(traj)
      sim = sims[0]
      for m in range(sim.num):
         traj = sim.get(m)
         values = sim.extract(traj)
         for i in range(values.shape[1]):
            mpl.subplot(values.shape[1], 1,i+1)
            mpl.plot(values[:,i])
            mpl.ylabel(sim.variables[i].name)
      mpl.show()
      self._finish_plot()


class Variance(Plot):
   def plot(self, sims, truth):
      traj = truth.get(0)
      values = truth.extract(traj)
      sim = sims[0]
      for m in range(sim.num):
         traj = sim.get(m)
         values = sim.extract(traj)
         for i in range(values.shape[1]):
            mpl.subplot(values.shape[1], 1,i+1)
            mpl.plot(values[:,i])
            mpl.ylabel(sim.variables[i].name)
      mpl.show()
      self._finish_plot()
