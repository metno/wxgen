from __future__ import division
import numpy as np
import datetime
import wxgen.util


class ClimateModel(object):
    """
    Base class for representing a climate state that can be used as external forcing for the weather
    generator.
    """
    def get(self, unixtimes):
        """ Returns a representation of the state for a given date

        Arguments:
           unixtimes (np.array): An array of unix times

        Returns:
           np.array: A 2D array representing the state. The first dimension equals the length of
              unixtimes, the second the number of variables in the state. Must be 2D even if the
              second second dimension is only 1.
        """
        raise NotImplementedError()


class Zero(ClimateModel):
    """ No climate forcing """
    def __init__(self):
        pass

    def get(self, unixtimes):
        return np.array([0 for t in unixtimes])


class Bin(ClimateModel):
    """
    State is determined by which bin the day of the year falls in. For a bin size of 10, then Jan
    1-10 are in bin 0, Jan 11-20 are in bin 1, and so forth.
    """
    def __init__(self, num_days):
        """
        Arguments;
           num_days (int): Number of days in each bin
        """
        self._num_days = num_days
        self._min_bin_fraction_size = 0.5

        # Ensure that you don't get a tiny bin at the end of the year. Instead, pool the last few days
        # together with the second last bin. For example, when using -b 30, the last 6 days of the
        # year will be in a bin of its own. Prevent this, by only creating a separate last bin if the
        # last bin is at least half the size of the other bins.
        max_bin = np.floor(365 / self._num_days)
        last_bin_start = self._num_days * max_bin
        last_bin_size = 365 - last_bin_start + 1
        min_bin_size = np.floor(self._num_days * self._min_bin_fraction_size)
        self._remove_last_bin = last_bin_size < min_bin_size

    def get(self, unixtimes):
        # Find the day of year on the interval [0, 365]
        day = np.array([wxgen.util.day_of_year(unixtime)-1 for unixtime in unixtimes])
        bin = day // self._num_days

        if self._remove_last_bin:
            # Find all times that are in the highest bin and put it in the previous bin
            max_bin = np.floor(365 / self._num_days)
            last_bin_start = self._num_days * max_bin
            I = day >= last_bin_start
            if max_bin == 0:
                # This should never happen, but lets just be safe just in case
                bin[I] = 0
            else:
                bin[I] = max_bin - 1

        bin = np.expand_dims(bin, 1)
        return bin


class Index(ClimateModel):
    """
    Use a climate index from a file
    """
    def __init__(self, filename, num_days=365):
        self._filename = filename
        self._num_days = num_days
        self._edges = np.array([-100, -1, 1, 100])
        fid = open(self._filename, 'r')
        self._index = dict()
        prev = 0
        for line in fid:
            words = line.strip().split(' ')
            words = [word for word in words if word != ""]
            year = int(words[0])
            month = int(words[1])
            day = int(words[2])
            date = int("%04d%02d%02d" % (year, month, day))
            unixtime = wxgen.util.date_to_unixtime(date)
            if words[3] == "NA":
                value = prev
            else:
                value = float(words[3])
                prev = value
            I = np.where(value > self._edges)[0]
            self._index[unixtime] = I[-1]

    def get(self, unixtimes):
        day = np.array([wxgen.util.day_of_year(unixtime)/self._num_days for unixtime in unixtimes])
        index = np.nan*np.zeros(len(unixtimes))
        for i in range(0, len(unixtimes)):
            unixtime = unixtimes[i]
            if unixtime in self._index:
                index[i] = self._index[unixtime]
            else:
                print("Missing: %d" % unixtime)
        bin = day + index * 100
        bin = np.expand_dims(bin, 1)
        return bin


class Combo(ClimateModel):
    """
    Combines several climate models. A matching state is one that matches the states from all models
    """
    def __init__(self, models):
        self.models = models

    def get(self, unixtimes):
        states = self.models[0].get(unixtimes)
        for m in range(1, len(self.models)):
            states = np.append(states, self.models[m].get(unixtimes), axis=1)

        return states
