import numpy as np


class Trajectory(object):
    """
    Represents a trajectory of states. This is represented as a sequence of indicies into the
    database: (Segment index, lead-time index).

    Attributes:
       indices (np.array): Indices for database
       length (int): Length of trajectory
    """
    def __init__(self, indices):
        """
        Arguments:
           indices (np.array): Array with dimensions [Segment index, lead-time index]
        """
        assert(len(indices.shape) == 2)
        assert(indices.shape[1] == 2)
        self.indices = indices

    @property
    def length(self) -> int:
        return len(self)

    def __str__(self):
        max_print = 10
        idx_str = ','.join(["[%d,%d]" % (index[0], index[1]) for index in self.indices[:max_print]])
        if len(self) <= max_print:
            return f"Trajectory(length={len(self)}, indices=[{idx_str}])"
        else:
            return f"Trajectory(length={len(self)}, indices=[{idx_str}, ...])"

    def __len__(self) -> int:
        return self.indices.shape[0]
