import inspect
import netCDF4
import matplotlib.pyplot as mpl
import numpy as np
import sys
import wxgen.util
import wxgen.transform
import wxgen.aggregator
import matplotlib.dates
import scipy.ndimage
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
   """
   Class for representing a verification plot of trajectory information
   """
   def __init__(self):
      self.filename = None
      self.dpi = 100
      self.fig_size = [10, 5]
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
      self.line_colors = None
      self.line_styles = None
      self.default_colors = ['r', 'b', 'g', [1, 0.73, 0.2], 'k']
      self.default_lines = ['-', '-', '-', '--']
      self.default_markers = ['o', '', '.', '']
      self.transform = wxgen.transform.Nothing()
      self.aggregator = wxgen.aggregator.Mean()
      self.clim = None
      self.cmap = None
      self.lat = None
      self.lon = None
      self.timescale = 1

   def plot(self, sims, truth):
      """
      Arguments:
         sims (list): List of simulation databases
         truth (wxgen.database): Truth database
      """
      raise NotImplementedError()

   def create_yearly_series(self, array):
      # Create 1-year long segments
      N = int(np.ceil(len(array)/365))
      array2 = np.zeros([365, N])
      for i in range(0, N):
         I = range(i*365, (i+1)*365)
         array2[:, i] = array[I]
      return array2

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
               xticklabels = [t.get_text() for t in ax.get_xticklabels()]
            ax.set_xscale("log")
            if self._sets_xticks:
               ax.set_xticks(xticks)
               ax.set_xticklabels(xticklabels)
         if self.ylog:
            if self._sets_yticks:
               yticks = ax.get_yticks()
               yticklabels = [t.get_text() for t in ax.get_yticklabels()]
            ax.set_yscale("log")
            if self._sets_yticks:
               ax.set_yticks(yticks)
               ax.set_yticklabels(yticklabels)
      if self.filename is None:
         mpl.show()
      else:
         mpl.gcf().set_size_inches(self.fig_size[0], self.fig_size[1])
         mpl.savefig(self.filename, bbox_inches='tight', dpi=self.dpi)

   def _get_color(self, i, total):
      """
      Returns a color specification (e.g. 0.3,0.3,1) that can be used in mpl to specify line color.
      Determined by looping through a database (self.line_colors). Returns the color for the i'th
      line in a plot of 'total' number of lines.

      _get_color together with _get_style can be used to specify unique color/style combinations for
      many lines. Color is cycled first, then style. I.e. the following order is default:
      r-o, b-o, g-o, ..., r-, b-, g-, ...

      Arguments:
         i (int): Which line is this?
         total (int): Total number of lines in plot

      Returns:
         list: A list of RGB values
      """
      if self.line_colors is not None:
         firstList = self.line_colors.split(",")
         numList = []
         finalList = []

         for string in firstList:
            if "[" in string:   # for rgba args
               if not numList:
                  string = string.replace("[", "")
                  numList.append(float(string))
               else:
                  wxgen.util.error("Invalid rgba arg \"{}\"".format(string))

            elif "]" in string:
               if numList:
                  string = string.replace("]", "")
                  numList.append(float(string))
                  finalList.append(numList)
                  numList = []
               else:
                  wxgen.util.error("Invalid rgba arg \"{}\"".format(string))

            # append to rgba lists if present, otherwise grayscale intensity
            elif wxgen.util.is_number(string):
               if numList:
                  numList.append(float(string))
               else:
                  finalList.append(string)

            else:
               if not numList:  # string args and hexcodes
                  finalList.append(string)
               else:
                  wxgen.util.error("Cannot read color args.")
         self.colors = finalList
         return self.colors[i % len(self.colors)]

      # use default colours if no colour input given
      else:
         self.colors = self.default_colors
         return self.colors[i % len(self.default_colors)]

   def _get_style(self, i, total, connectingLine=True, lineOnly=False):
      """ Returns a string (e.g. -o) that can be used in mpl to specify line
      style. Determined by looping through a database (self.line_styles).
      Returns the style for the i'th line in a plot of 'total' number of lines.

      Arguments:
         i (int): Which line is this?
         total (int): Total number of lines in plot
         connectingLine: If True, add a connecting line (e.g. -o) between the
            markers.  Otherwise only a marker will be used (e.g. o)
         lineOnly: If True, don't include the marker (e.g. -)
      """
      if self.line_styles is not None:
         listStyles = self.line_styles.split(",")
         # loop through input linestyles (independent of colors)
         I = i % len(listStyles)
         return listStyles[I]

      else:  # default linestyles
         I = (i / len(self.colors)) % len(self.default_lines)
         line = self.default_lines[I]
         marker = self.default_markers[I]
         if lineOnly:
            return line
         if connectingLine:
            return line + marker
         return marker

   def _plot_truth(self, x, y, isCont=True, zorder=0, label="Truth"):
      if isCont:
         mpl.plot(x, y, ".-", color="gray", lw=5, label=label, zorder=zorder)
      else:
         mpl.plot(x, y, "o", color="gray", ms=self.ms, label=label,
               zorder=zorder)


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
                  mpl.subplot(X, Y, index)
                  Ivar = Ivars[i]
                  mpl.plot(values[:, Ivar], '-')
                  mpl.ylabel(sim.variables[Ivar].name)
                  mpl.title(sim.name)
                  mpl.xlabel("Time (days)")
                  mpl.grid()

      if truth is not None:
         traj = truth.get(0)
         values = truth.extract(traj)
         for i in range(len(Ivars)):
            index = (X-1)*Y+i+1
            mpl.subplot(X, Y, index)
            Ivar = Ivars[i]
            mpl.plot(values[:, Ivar], lw=5, color="red")
            mpl.ylabel(truth.variables[Ivar].name)
            mpl.grid()

      self._finish_plot()


class Histogram(Plot):
   """ Plot histogram for annual values for a given statistic """
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

      for i in range(len(Ivars)):
         Ivar = Ivars[i]
         xlim = None
         if sims is not None:
            for s in range(len(sims)):
               index = s*Y+i+1
               mpl.subplot(X, Y, index)
               sim = sims[s]
               agg = np.zeros(0)
               for m in range(sim.num):
                  traj = sim.get(m)
                  values = sim.extract(traj)[:, Ivar]
                  values = self.transform(values)
                  if len(values) < 365:
                     wxgen.util.error("Simulations must be longer than 365 days long")
                  values = self.create_yearly_series(values)
                  curr_agg = np.zeros(values.shape[1])
                  for k in range(values.shape[1]):
                     curr_agg[k] = self.aggregator(values[:, k])
                  agg = np.append(agg, curr_agg)

               if self.thresholds is not None:
                  mpl.hist(agg, self.thresholds, normed=1)
               else:
                  mpl.hist(agg, normed=1)
               mpl.ylabel("Density")
               mpl.xlabel(sim.variables[Ivar].name)
               mpl.title(sim.name)
               if self.thresholds is not None:
                  dx = self.thresholds[1]-self.thresholds[0]
                  mpl.ylim([0, 1.0/dx])

         if truth is not None:
            index = (X-1)*Y+i+1
            mpl.subplot(X, Y, index)
            traj = truth.get(0)
            values = truth.extract(traj)[:, Ivar]
            values = self.transform(values)   # 5 years
            values = self.create_yearly_series(values)  # 5 x 1 year
            agg = np.zeros(values.shape[1])
            for k in range(values.shape[1]):
               agg[k] = self.aggregator(values[:, k])
            if self.thresholds is not None:
               mpl.hist(agg, self.thresholds, normed=1)
            else:
               mpl.hist(agg, normed=1)
            mpl.ylabel("Density")
            mpl.xlabel(truth.variables[Ivar].name)
            mpl.title(truth.name)
            if self.thresholds is not None:
               dx = self.thresholds[1]-self.thresholds[0]
               mpl.ylim([0, 1.0/dx])

      self._finish_plot()


class Variance(Plot):
   def __init__(self):
      Plot. __init__(self)
      self._sets_xticks = True
      self._normalize = False
      self._normalization_window = 11

   def plot(self, sims, truth):
      if self.thresholds is None:
         scales = [1, 3, 7, 11, 31, 61, 181, 365]
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
            values = truth.extract(traj)[:, Ivar]
            # I = np.where(np.isnan(values) == 0)[0]
            # values = values[I]
            truth_var = self.compute_truth_variance(values, scales)
            self._plot_truth(scales, truth_var, label="Truth")
            # mpl.plot(scales, truth_var, 'o-', label="Truth")
            mpl.title(truth.variables[Ivar].name)
            if self._normalize:
               mpl.ylabel("Normalized variance")
            else:
               mpl.ylabel("Variance ($%s^s$)" % truth.variables[Ivar].units)

         if sims is not None:
            for s in range(len(sims)):
               sim = sims[s]
               sim_values = np.zeros([sim.length, sim.num])
               col = self._get_color(s, len(sims))
               for m in range(sim.num):
                  traj = sim.get(m)
                  q = sim.extract(traj)
                  sim_values[:, m] = q[:, Ivar]
               sim_var = self.compute_sim_variance(sim_values, scales)
               mpl.plot(scales, sim_var, 'o-', label=sim.name, color=col)
         ticks = np.array([1, 7, 30, 365])
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
      import astropy.convolution
      # Create 1-year long segments
      truth = self.create_yearly_series(array)

      """
      Remove climatology so we can look at annomalies. Use separate obs and fcst climatology
      otherwise the fcst variance is higher because obs gets the advantage of using its own
      climatology.
      """
      clim = wxgen.util.climatology(truth, self._normalization_window)
      std = 1
      if self._normalize:
         std = np.nanstd(truth, axis=1)

      N = truth.shape[1]
      for i in range(0, N):
         truth[:, i] = (truth[:, i] - clim)/std

      # Compute variance
      variance = np.zeros(len(scales))
      for i in range(0, len(scales)):
         s = scales[i]
         c = [1.0/s] * s
         truth_c = np.zeros([truth.shape[0], N], float)
         for e in range(0, N):
            truth_c[:, e] = astropy.convolution.convolve(truth[:, e], 1.0/s*np.ones(s))
         if s > 1:
            truth_c = truth_c[(s/2):(-s/2+1), :]
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
      import astropy.convolution
      N = array.shape[1]

      clim = wxgen.util.climatology(array, self._normalization_window)
      std = 1
      if self._normalize:
         std = np.nanstd(array, axis=1)
      values = copy.deepcopy(array)
      for i in range(0, N):
         values[:, i] = (values[:, i] - clim)/std

      # Compute variance
      variance = np.nan*np.zeros(len(scales))
      for i in range(0, len(scales)):
         s = scales[i]
         if array.shape[0] >= s:
            c = [1.0/s] * s
            sim_c = np.zeros([values.shape[0], N], float)
            for e in range(0, N):
               sim_c[:, e] = astropy.convolution.convolve(values[:, e], 1.0/s*np.ones(s))
            if s > 1:
               sim_c = sim_c[(s/2):(-s/2+1), :]
            variance[i] = np.nanvar(sim_c)
      return variance


class Map(Plot):
   def plot(self, sims, truth):
      if truth is not None:
         sims += [truth]
      if self.vars is None:
         Ivars = range(len(sims[0].variables))
      else:
         Ivars = self.vars

      X = len(sims)
      Y = len(Ivars)
      import mpl_toolkits.basemap
      for v in range(len(Ivars)):
         Ivar = Ivars[v]
         for s in range(len(sims)):
            count = 0
            index = s*Y+v+1
            mpl.subplot(X, Y, index)
            sim = sims[s]
            sim_values = np.zeros([sim.Y, sim.X])
            lats = sim.lats
            lons = sim.lons
            dlat = 1.0
            dlon = 1.0
            llcrnrlat = max(-90, np.min(lats) - dlat / 10)
            urcrnrlat = min(90, np.max(lats) + dlat / 10)
            llcrnrlon = np.min(lons) - dlon / 10
            urcrnrlon = np.max(lons) + dlon / 10
            res = "l"
            map = mpl_toolkits.basemap.Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                  urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, projection='cyl',
                  resolution=res)
            for m in range(sim.num):
               traj = sim.get(m)
               q = sim.extract_grid(traj)
               sim_values += self.aggregator(self.transform(q[:, :, :, Ivar]), axis=0)
               count += 1
            [x, y] = map(lons, lats)
            if self.clim is not None:
               map.contourf(x, y, sim_values/count, np.linspace(self.clim[0], self.clim[1], 11), label=sim.name, cmap=self.cmap)
            else:
               map.contourf(x, y, sim_values/count, label=sim.name, cmap=self.cmap)
            cb = map.colorbar()
            if self.clim is not None:
               mpl.clim(self.clim)
               cb.set_clim(self.clim)
            map.drawcoastlines(linewidth=1)
            map.drawcountries(linewidth=2)
            map.drawmapboundary()
            dy = 1
            dx = 1
            map.drawparallels(np.arange(-90., 120., dy), labels=[1, 0, 0, 0])
            map.drawmeridians(np.arange(-180., 420., dx), labels=[0, 0, 0, 1])
            map.fillcontinents(color=[0.7, 0.7, 0.7], zorder=-1)
            mpl.title(sim.name)
      self._finish_plot()


class Jump(Plot):
   def plot(self, sims, truth):
      if truth is not None:
         sims += [truth]
      if self.vars is None:
         Ivars = range(len(sims[0].variables))
      else:
         Ivars = self.vars

      X = 1
      Y = len(Ivars)
      db_length = 10
      for v in range(len(Ivars)):
         Ivar = Ivars[v]
         index = v+1
         mpl.subplot(X, Y, index)
         for s in range(len(sims)):
            count = 0
            sim = sims[s]
            values = np.zeros([db_length-1])
            counts = np.zeros([db_length-1])
            for m in range(sim.num):
               traj = sim.get(m)
               q = sim.extract_grid(traj)
               for i in range(0, sim.length-1):
                  I = i % (db_length-1)
                  curr = np.mean(np.abs(q[i, :, :, Ivar] - q[i+1, :, :, Ivar]))
                  values[I] += curr
                  counts[I] += 1
            values = values / counts
            col = self._get_color(s, len(sims))
            mpl.plot(np.arange(0.5, db_length - 0.5), values, '-o', color=col, label=sim.name)
         mpl.xlabel("Lead time (days)")
         mpl.ylabel("Average absolute jump")
         mpl.legend(loc="best")
         mpl.grid()
         mpl.ylim(ymin=0)
      self._finish_plot()


class TimeStat(Plot):
   def plot(self, sims, truth):
      if truth is not None:
         sims += [truth]
      if self.vars is None:
         Ivars = range(len(sims[0].variables))
      else:
         Ivars = self.vars

      X = 1
      Y = len(Ivars)
      for v in range(len(Ivars)):
         Ivar = Ivars[v]
         index = v+1
         mpl.subplot(X, Y, index)
         for s in range(len(sims)):
            count = 0
            sim = sims[s]
            if self.timemod is None:
               L = sim.length
            else:
               L = self.timemod

            """
            Take values from all gridpoints, leadtimes, and members and then aggregate at the end.
            This could potentially require a large memory footprint, but since the aggregation isn't
            always just a sum, we have to load all data first, then aggregate (instead of
            iteratively summing up values, without keeping all values in memory at once),
            """
            values = [np.zeros([0])]*L
            for m in range(sim.num):
               traj = sim.get(m)
               q = sim.extract_grid(traj)
               for i in range(L):
                  I = range(i, q.shape[0], L)
                  values[i] = np.append(values[i], self.transform(q[I, :, :, Ivar]).flatten())
            values_agg = np.zeros(L)
            for i in range(L):
               values_agg[i] = self.aggregator(values[i])
            col = self._get_color(s, len(sims))
            mpl.plot(range(L), values_agg, '-o', color=col, label=sim.name)
         mpl.xlabel("Lead time (days)")
         mpl.ylabel("%s" % (self.aggregator.name().capitalize()))
         mpl.legend(loc="best")
         mpl.grid()
      self._finish_plot()


class SortStat(Plot):
   def plot(self, sims, truth):
      if truth is not None:
         sims += [truth]
      if self.vars is None:
         Ivars = range(len(sims[0].variables))
      else:
         Ivars = self.vars

      X = len(sims)
      Y = len(Ivars)
      for v in range(len(Ivars)):
         Ivar = Ivars[v]
         for s in range(len(sims)):
            index = v+1
            index = s*Y+v+1
            mpl.subplot(X, Y, index)
            count = 0
            sim = sims[s]
            if self.timemod is None:
               L = sim.length
            else:
               L = self.timemod
            if L > 100:
               wxgen.util.error("Too many lines. Consider -tm.")
            values = [np.zeros([0])]*L
            for m in range(sim.num):
               traj = sim.get(m)
               q = sim.extract_grid(traj)
               for i in range(L):
                  I = range(i, q.shape[0], L)
                  values[i] = np.append(values[i], self.transform(q[I, :, :, Ivar]).flatten())
            for i in range(L):
               col = self._get_color(i, L)
               mpl.plot(np.sort(values[i]), np.linspace(0, 1, len(values[i])), '-o', color=col, label="Day %d" % i, ms=0)
            mpl.xlabel(sim.variables[Ivar].name)
            mpl.ylabel("Quantile")
            mpl.title(sim.name)
            mpl.legend(loc="best")
            mpl.grid()
      self._finish_plot()


class CovarMap(Plot):
   def plot(self, sims, truth):
      if self.lat is None or self.lon is None:
         wxgen.util.error("-lat and/or -lon not specified")

      if truth is not None:
         sims += [truth]
      if self.vars is None:
         Ivars = range(len(sims[0].variables))
      else:
         Ivars = self.vars

      Xref = 5
      Yref = 10
      Y = len(sims)
      X = len(Ivars)
      import mpl_toolkits.basemap
      import astropy.convolution
      for v in range(len(Ivars)):
         Ivar = Ivars[v]
         for s in range(len(sims)):
            count = 0
            index = s*X+v+1
            mpl.subplot(X, Y, index)
            sim = sims[s]
            sim_values = np.zeros([sim.Y, sim.X])
            lats = sim.lats
            lons = sim.lons
            dlat = 1.0
            dlon = 1.0
            llcrnrlat = max(-90, np.min(lats) - dlat / 10)
            urcrnrlat = min(90, np.max(lats) + dlat / 10)
            llcrnrlon = np.min(lons) - dlon / 10
            urcrnrlon = np.max(lons) + dlon / 10
            [Xref, Yref] = wxgen.util.get_i_j(lats, lons, self.lat, self.lon)
            res = "l"
            map = mpl_toolkits.basemap.Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                  urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, projection='cyl',
                  resolution=res)
            for m in range(sim.num):
               traj = sim.get(m)
               q = sim.extract_grid(traj)
               val = q[:, :, :, Ivar]
               ref = np.swapaxes(np.tile(q[:, Xref, Yref, Ivar], [q.shape[2], q.shape[1], 1]), 0, 2)
               scale = self.timescale
               if scale % 2 != 1:
                  wxgen.util.error("Time scale must be an odd number")
               if scale > 1:
                  sarray = 1.0/scale*np.ones(scale)
                  for i in range(0, ref.shape[1]):
                     for j in range(0, ref.shape[2]):
                        ref[:, i, j] = astropy.convolution.convolve(ref[:, i, j], sarray)
                        val[:, i, j] = astropy.convolution.convolve(val[:, i, j], sarray)

                  # Remove edges in convolution
                  ref = ref[(scale/2):(-scale/2+1), :, :]
                  val = val[(scale/2):(-scale/2+1), :, :]

               sim_values += wxgen.util.correlation(val, ref, axis=0)
               count += 1
            [x, y] = map(lons, lats)
            if self.clim is not None:
               map.contourf(x, y, sim_values/count, np.linspace(self.clim[0], self.clim[1], 11),
                     label=sim.name, cmap=self.cmap, extend="both")
            else:
               map.contourf(x, y, sim_values/count, label=sim.name, cmap=self.cmap, extend="both")
            mpl.plot(x[Xref, Yref], y[Xref, Yref], '*', ms=15, mfc="yellow", mec="k")
            cb = map.colorbar()
            if self.clim is not None:
               mpl.clim(self.clim)
               cb.set_clim(self.clim)
            map.drawcoastlines(linewidth=1)
            map.drawcountries(linewidth=2)
            map.drawmapboundary()
            dy = 1
            dx = 1
            map.drawparallels(np.arange(-90., 120., dy), labels=[1, 0, 0, 0])
            map.drawmeridians(np.arange(-180., 420., dx), labels=[0, 0, 0, 1])
            map.fillcontinents(color=[0.7, 0.7, 0.7], zorder=-1)
            mpl.title(sim.name)
      self._finish_plot()
