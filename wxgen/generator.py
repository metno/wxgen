from typing import Optional
import numpy as np
from wxgen.database import Database, SegmentIndices
import wxgen.metric
from wxgen.trajectory import Trajectory
import wxgen.util
import wxgen.climate_model
from tqdm import tqdm

import logging

logger = logging.getLogger(__name__)

class Generator(object):
    """
    This class generates long trajectories from segments in a database
    """

    def __init__(self, database: Database, metric=wxgen.metric.Rmsd(), model=None):
        self._database = database
        self._metric = metric
        self.prejoin = None
        self.policy = "random"
        self.stagger = False
        self.start_unixtime = wxgen.util.date_to_unixtime(20170101)
        self.db_start_date = None
        self.db_end_date = None

    def get(self, N: int, T: int, initial_state: Optional[np.ndarray] = None) -> list[Trajectory]:
        """Randomly generated trajectories

        If initial_state is provided then the trajectory will start with a state similar to this. If
        no initial state is provided, start with a random segment from the database.

        Arguments:
           N: Number of trajectories
           T: Number of timesteps in each trajectory
           initial_state: An array of the initial state (must be of length V)

        Returns:
            List of N trajectories, where each trajectory has a length of T
        """
        trajectories = list()
        if self._database.length < 2:
            raise RuntimeError("Cannot create simulation with a database with segments shorter than 2 timesteps")

        for n in range(0, N):
            logger.info("Generating trajectory %d/%d", n+1, N)
            trajectory_indices = -1 + np.zeros([T, 2], int)

            time = self.start_unixtime
            time_of_day = self.start_unixtime % 86400

            if initial_state is None:
                state_curr = None
            else:
                state_curr = initial_state

            """
            Assemble a trajectory by concatenating appropriate segments. Start by finding a segment
            that has a starting state that is similar to the requested initial state. When repeating,
            overwrite the end state of the previous segment. This means that if the segment is 10 days
            long, we are only using 9 days of the segment.

            If the last segment fits exactly, we do not need to search for a new state to fill in the
            last day with new data. This is the case when start == T-1.
            """
            start = 0  # Starting index into output trajectory where we are inserting a segment
            join = 0
            with tqdm(total=T, mininterval=1) as pbar:
                while start < T-1:
                    inc = start - pbar.n
                    pbar.update(inc)

                    """
                    Prejoin multiple segments that are nearby in time. This is done by passing
                    'search_times' to get_random.
                    """
                    search_times = None
                    if join > 0:
                        end_times = self._database.inittimes[segment_curr.indices[-1, 0]] + segment_curr.indices[-1, 1]*self._database.timestep
                        search_times = [end_times - 5*86400, end_times + 5*86400]
                    logger.debug("Found random segment")
                    if state_curr is None:
                        logger.debug("Target state: None")
                    else:
                        logger.debug("Target state: %s" % ' '.join(["%0.2f" % x for x in state_curr]))
                    state_curr_unix_time = time
                    segment_curr = self.get_random(state_curr, time_of_day, self._metric, state_curr_unix_time, search_times)
                    timesteps_per_day = int(86400 / self._database.timestep)

                    """
                    Stagger the trajectories so that they don't all potentially have jumps at the same
                    leadtimes. This is done by truncating the first segment to a random length. Note that
                    the upper end of randint is exclusive, hence the "+ 1".

                    For subdaily timesteps, this must cut whole days
                    """
                    if self.stagger and start == 0:
                        I = np.random.randint(1, segment_curr.indices.shape[0]/timesteps_per_day + 1) * timesteps_per_day
                        if I > segment_curr.indices.shape[0] + 1:
                            I = segment_curr.indices.shape[0] + 1
                        segment_curr.indices = segment_curr.indices[0:I, :]

                    indices_curr = segment_curr.indices

                    """
                    Account for the fact that the desired trajectory length is not a whole multiple of the
                    segment length: Only take the first part of the segment if needed.

                    Also account for the fact that the last timestep in the segment must be at the same time
                    of day as the first timestep in the segment so that matching occurrs with the same time
                    of day (and overwrite the necessary time-step in question).
                    """
                    Tsegment = len(indices_curr)
                    end = start + Tsegment  # Ending index
                    end = min(end, T)  # If this is the last segment, then make sure it doesn't go past the length of the desired trajectory
                    if start == 0:
                        end = end-1
                        Iout = range(start, end)  # Index into trajectory
                        Iin = range(0, end - start)  # Index into segment
                    else:
                        Iout = range(start, end)  # Index into trajectory
                        Iin = range(0, end - start)  # Index into segment
                    trajectory_indices[Iout, :] = indices_curr[Iin, :]
                    # print Tsegment, start, end, time_of_day//3600
                    # print Iin, Iout

                    logger.debug("Current state: %s", state_curr)
                    logger.debug("Chosen segment: %s", segment_curr)
                    logger.debug("Trajectory indices: %s", Iout)
                    logger.debug("Segment indices: %s", Iin)

                    # Get the last state of the segment, and corresponding time
                    state_curr = self._database.extract_matching(segment_curr)[Iin[-1], :]
                    start = start + Tsegment - 1
                    time = time + (Tsegment - 1) * self._database.timestep
                    time_of_day = time % 86400
                    if self.prejoin is not None and self.prejoin > 0:
                        join = (join + 1) % self.prejoin

                if len(np.where(trajectory_indices == -1)[0]) > 0:
                    raise RuntimeError("Internal error. The trajectory was not properly filled")
                trajectory = wxgen.trajectory.Trajectory(trajectory_indices)
                if logger.debug:
                    logger.debug("Trajectory: %s", trajectory)
                trajectories.append(trajectory)

            return trajectories

    def get_random(self, target_state: Optional[np.ndarray], 
                   time_of_day: int, 
                   metric: wxgen.metric.Metric, 
                   target_state_unix_time: Optional[int], 
                   time_range: Optional[list] = None) -> Trajectory:
        """
        Returns a pseudo-random segment from the database chosen based on weights computed by a metric

        Arguments:
           target_state: Try to match this state when finding the trajectory. One value
              for each variable in the database. Or if None, then pick a random start segment.
           time_of_day: Time of day to start simulation
           metric: Metric to use when finding matches
           target_state_unix_time: Time stamp of the target state.
           time_range: Start and end unixtimes for the search. If None, then do not restrict.

        Returns:
           Random trajectory based on weights computed by a metric
        """
        assert(target_state is None or np.sum(np.isnan(target_state)) == 0)
        assert(self._database._data_matching.shape[2] == self._database.num)

        # TODO: move idx0 into cutting the weights / segment indices?
        idx0 = np.searchsorted(self._database.leadtimes, time_of_day)

        weights = self._weights_per_state(target_state=target_state, idx0=idx0, metric=metric)

        # Find valid segments
        idx_segments, do_prejoin = self._find_valid_segments(weights, time_range)
        # TODO: why shouldn't we do this when prejoining? 
        if (not do_prejoin):
            assert target_state_unix_time is not None, f"{target_state_unix_time=} has to be provided"
            idx_segments = self._filter_on_climate_state(target_state_unix_time=target_state_unix_time, 
                                                         idx_segments=idx_segments)

        weights = self._flip_weights_for_negative_metric(weights, metric)
        tr = self._get_trajectory(weights=weights, idx_segments=idx_segments, idx0=idx0)
        return tr

    def idx_not_nan_1d(self, arr: np.ndarray) -> np.ndarray:
        """Index of not-nan elements in arr"""
        assert arr.ndim == 1
        return (~np.isnan(arr)).nonzero()[0]

    def idx_where_true_1d(self, arr: np.ndarray) -> np.ndarray:
        """Index of elements that are `True` in arr"""
        assert arr.ndim == 1
        return arr.nonzero()[0]
    
    def _find_valid_segments(self, weights: np.ndarray, time_range: Optional[list]) -> tuple[SegmentIndices, bool]:
        """Find valid segements, so weight is not nan, and that are within the specified time range.
        Args:
            weights: weight for each segment in the database
            time_range: Start and end unixtimes for the search. If None, then do not restrict.

        Returns: 
            Returns tuple: (Array of indices of valid segments, do_prejoin)
        """
        do_prejoin = False
        if time_range is None and self.db_start_date is None and self.db_end_date is None:
            idx_segments = self.idx_not_nan_1d(weights)
        elif self.db_start_date is not None and self.db_end_date is not None:
            db_start_date = wxgen.util.date_to_unixtime(self.db_start_date)
            db_end_date = wxgen.util.date_to_unixtime(self.db_end_date)
            cond = (np.isnan(weights) == 0) & (self._database.inittimes > db_start_date) & (self._database.inittimes < db_end_date)
            idx_segments = self.idx_where_true_1d(cond)
        elif self.db_start_date is not None:
            db_start_date = wxgen.util.date_to_unixtime(self.db_start_date)
            cond = (np.isnan(weights) == 0) & (self._database.inittimes > db_start_date)
            idx_segments = self.idx_where_true_1d(cond)
        elif self.db_end_date is not None:
            db_end_date = wxgen.util.date_to_unixtime(self.db_end_date)
            cond = (np.isnan(weights) == 0) & (self._database.inittimes < db_end_date)
            idx_segments = self.idx_where_true_1d(cond)
        else:
            do_prejoin = True
            cond = (np.isnan(weights) == 0) & (self._database.inittimes > time_range[0]) & (self._database.inittimes < time_range[1])
            idx_segments = self.idx_where_true_1d(cond)
            if len(idx_segments) == 0:
                date_range = [wxgen.util.unixtime_to_date(t) for t in time_range]
                logger.warning("Skipping this prejoin: No valid segment that start in date range [%d, %d]",
                      date_range[0], date_range[1])
                idx_segments = self.idx_not_nan_1d(weights)
                # Without any prejoin segments, revert to the original plan of just finding a random segment
                do_prejoin = False
        return idx_segments, do_prejoin
    
    def _weights_per_state(self, target_state: np.ndarray, idx0: int, metric: wxgen.metric.Metric) -> np.ndarray:
        """Weights for each segment in the database chosen based on weights computed by a metric

        Args:
            target_state: Try to match this state when finding the trajectory
            idx0: First index of the database to take into account 
            metric: Metric to use when finding matches

        Returns:
            Weight for each segment in the database

        """
        # TODO: should climat state matching be optional?
        # TODO: climate state can also mean here that one groups into bins of +-x days, roughyl

        if target_state is None:
            """
            A bit of a hack. Here we want to select a random segment, therefore set all the
            weights to one.  However, we do not want to select a segment with missing values, so
            we have to insert nans for such segments.
            """
            weights = np.ones(self._database.num)
            idx_segment_is_nan = np.sum(np.isnan(self._database._data_matching[idx0, :, :]), axis=0) > 0
            weights[idx_segment_is_nan] = np.nan
        else:
            weights = metric.compute(target_state, self._database._data_matching[idx0, :, :])
        return weights
    
    def _filter_on_climate_state(self, target_state_unix_time: int,
                                idx_segments: SegmentIndices) -> SegmentIndices:
        """Filter on 'climate state' (if applies)
        Args:
            target_state_unix_time: Unix time stamp of the 'Climate state' that should be matched
            idx_segments: candidates segment member indices (before filter)
        
        Returns:
            Filtered segment member indices

        """
        logger.debug("Filter on climate states")
        time_stamps_candidates = self._database.inittimes[idx_segments]
        idx_filtered = self._database.model.match(target_state_unix_time, time_stamps_candidates)
        if len(idx_filtered) == 0:
            raise RuntimeError(f"Cannot find a segment with matching climate states for {target_state_unix_time}")
        idx_segments = idx_segments[idx_filtered]
        return idx_filtered
    
    def _flip_weights_for_negative_metric(self, weights: np.ndarray, metric: wxgen.metric.Metric) -> np.ndarray:
        """Flip weights if the metric is negative oriented"""
        if metric._orientation == -1:
            I0 = np.where(weights < 1e-3)[0]
            I1 = np.where(weights >= 1e-3)[0]
            # Ensure we do not get too high weights
            weights[I1] = 1.0/weights[I1]
            weights[I0] = 1e3
        return weights
    
    def _get_trajectory(self, weights: np.ndarray, idx_segments: SegmentIndices, idx0: int) -> Trajectory:
        """Random trajectory, chosen accoridng to the weights
        
        Args:
            weights: weight for each segment in the database
            idx_segments: segment member indices to consider
            idx0: First index of the database to take into account 

        Returns:
            Pseudo-randomly chosen trajectory
        """
        weights = weights[idx_segments]
        idx_in_segments_list = wxgen.util.random_weighted(weights, self.policy)
        idx_segment_selected = idx_segments[idx_in_segments_list]
        tr = self._database.get(idx_segment_selected, i_lead_time_start=idx0)

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Num candidates:  %d", len(weights))
            times_of_segements = self._database.inittimes[idx_segments]
            logger.debug("Date range:  %d %d", 
                         wxgen.util.unixtime_to_date(np.min(times_of_segements)), 
                         wxgen.util.unixtime_to_date(np.max(times_of_segements)))
            _tmp = self._database._data_matching[idx0, :, idx_segment_selected]
            logger.debug("Found state:  %s" % ' '.join(["%0.2f" % x for x in _tmp]))
            logger.debug("Found date: %s (%i)", 
                         wxgen.util.unixtime_to_date(self._database.inittimes[idx_segment_selected]), 
                         idx_segment_selected)
            logger.debug("Weight (max weight): %s (%s)", weights[idx_in_segments_list], np.max(weights))

        return tr