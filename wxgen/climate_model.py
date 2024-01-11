from __future__ import division
import numpy as np
import wxgen.util
from abc import ABC, abstractmethod
import shyft.time_series as sts

class ClimateModel(ABC):
    """
    Base class for representing a climate state that can be used as external forcing for the weather
    generator.
    """
    @abstractmethod
    def match(self, target: int, candidates: np.ndarray[int]) -> np.ndarray[int]:
        """Returns a representation of the state for a given date

        Args:
            target: Target states time stamp
            candidates: List of time stamps of possible matching state

        Returns:
           Indices of matching states
        """
        ...

class Zero(ClimateModel):
    """ No climate forcing """
    def __init__(self):
        pass

    def match(self, target: int, candidates: np.ndarray[int]) -> np.ndarray[int]:
        return np.arange(len(candidates))


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

    def _get(self, unixtimes):
        # Find the day of year on the interval [0, 365]
        day = wxgen.util.day_of_year(unixtimes)-1
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
    
    def match(self, target: int, candidates: np.ndarray[int]) -> np.ndarray[int]:
        target_state = self._get([target])[0]
        candidate_states = self._get(candidates)
        idx_match = np.where(target_state == candidate_states)[0]
        return idx_match


class CloseDayOfYear(ClimateModel):
    """
    The same climate state means having the same day of year +- a time-window. Similar to `Bin`, but 
    with a gliding time window. Does not take into account leap-years
    """
    def __init__(self, window_days):
        """
        Arguments;
           num_days (int): Number of days in the time windlow
        """
        if window_days > 180:
            raise ValueError("Using a time window of more than several month does not make sense.")
        self.window_days = window_days
        self._year_index = np.arange(365)

        self._candidates_doy_cache = None
        self._candidates_cache = None
        
    def match(self, target: int, candidates: np.ndarray[int]) -> np.ndarray[int]:
        t0 = sts.utctime_now()
        target_doy = wxgen.util.day_of_year([target])
        idx_range_selected = np.take(self._year_index, [target_doy - self.window_days, target_doy + self.window_days], mode="wrap")

        if (self._candidates_doy_cache is None) or (not np.array_equal(self._candidates_cache, candidates)):
            candidates_doy = wxgen.util.day_of_year(candidates)
            self._candidates_doy_cache = candidates_doy
            self._candidates_cache = candidates
        else:
            candidates_doy = self._candidates_doy_cache
    
        idx_first = idx_range_selected[0]
        idx_last = idx_range_selected[-1]
        if idx_first < idx_last:
            # "normal"
            is_match = (candidates_doy >= idx_first) & (candidates_doy <= idx_last)
        else:
            # wrap around occured, e.g. each for target_doy = 5, and time_window = 10 -> idxes = [360, 15]
            idx_first, idx_last = idx_last, idx_first
            is_match = (candidates_doy <= idx_first) | (candidates_doy >= idx_last)

        idx_match = np.where(is_match)[0]
        return idx_match

