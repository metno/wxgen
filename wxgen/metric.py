import inspect
import numpy as np
import sys
import wxgen.util


def get_all():
    """ Returns a list of all metric classes """
    temp = inspect.getmembers(sys.modules[__name__], inspect.isclass)
    return temp


def get(name):
    """ Returns a metric object of a class with the given name """
    metrics = get_all()
    m = None
    for mm in metrics:
        if(name == mm[0].lower()):
            m = mm[1]()
    if m is None:
        raise RuntimeError("Cannot find metric called '%s'" % name)
    return m


class Metric(object):
    """
    Class to compute the similarity between two states
    """
    # +1 means higher values are better, -1 means lower values are better
    _orientation = 1

    def __init__(self):
        pass

    def compute(self, state1, state2):
        """
        Compute the comparison score between the two states

        Arguments:
           state1 (np.array): An array of dimensions (V,) or (V,N)
           state2 (np.array): An array of dimensions (V,) or (V,N)

        If state1 is (V,N), then state2 must also be (V,N). If state1 is (V,) and state2 (V,N),
        then state1 ir repeated for all N.
        """
        if state1.shape != state2.shape:
            # Ensure states are the same size
            state1 = wxgen.util.resize(state1, state2.shape)

        return self._compute(state1, state2)

    def _compute(self, state1, state2):
        """
        state1 and state2 are guaranteed to be the same size
        """
        raise NotImplementedError()


class Rmsd(Metric):
    """ Root mean squared difference of the states

    Arguments:
       weights (np.array): array of variable-weights
    """
    _orientation = -1

    def __init__(self, weights=None):
        self._weights = weights
        if self._weights is None:
            self._weights = 0 # todo: remove this default?

        self._weights_reshaped_cache = {}

    def _weights_reshaped(self, repetitions: int) -> np.ndarray:
        if repetitions not in self._weights_reshaped_cache:
            arr = np.tile(self._weights, (repetitions,1)).T
            self._weights_reshaped_cache[repetitions] = arr
        return self._weights_reshaped_cache[repetitions]

    def _compute(self, state1, state2):
        # weights = wxgen.util.resize(0, state2.shape) # = old version
        weights = self._weights_reshaped(state2.shape[1])
        total = np.sum(weights*(state1 - state2)**2, axis=0)
        return np.sqrt(total)


class Max(Metric):
    """ Max difference of the states

    Arguments:
       weights (np.array): array of variable-weights
    """
    _orientation = -1

    def __init__(self, weights=None):
        self._weights = weights
        if self._weights is None:
            self._weights = 0

    def _compute(self, state1, state2):
        weights = wxgen.util.resize(self._weights, state2.shape)
        total = np.max(weights*abs(state1 - state2), axis=0)
        return total


class Mad(Metric):
    """ Mean absolute difference of the states

    Arguments:
       weights (np.array): array of variable-weights
    """
    _orientation = -1

    def __init__(self, weights=None):
        self._weights = weights
        if self._weights is None:
            self._weights = 0

    def _compute(self, state1, state2):
        weights = wxgen.util.resize(self._weights, state2.shape)
        total = np.mean(weights*abs(state1 - state2), axis=0)
        return total


class Exp(Metric):
    """ Exponential score based on exp(-sum(|factor * diff|))

    Arguments:
       factors (np.array): array of variable-weights
    """
    _orientation = 1

    def __init__(self, factors=None):
        self._factors = factors
        if factors is None:
            self._factors = 1

    def _compute(self, state1, state2):
        factors = wxgen.util.resize(self._factors, state2.shape)
        total = np.exp(-np.sum(factors*abs(state1 - state2), axis=0))
        return total
