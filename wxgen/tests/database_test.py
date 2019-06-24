import unittest
import wxgen.generator
import wxgen.database
import numpy as np


class DatabaseTest(unittest.TestCase):
    def test_simple(self):
        db = wxgen.database.Netcdf('wxgen/tests/files/test_simple.nc')
        truth = db.get_truth(start_date=20180101, end_date=20180102)
        np.testing.assert_array_equal(truth.indices, np.array([[0, 0], [0, 1], [0, 2], [0, 3], [0, 4]]))

        truth = db.get_truth(start_date=20180101, end_date=20180103)
        np.testing.assert_array_equal(truth.indices, np.array([[0, 0], [0, 1], [0, 2], [0, 3], [0, 4], [0, 5], [0, 6], [0, 7], [1, 0]]))

    def test_deacc(self):
        db = wxgen.database.Netcdf('wxgen/tests/files/test_simple.nc')
        db.deacc = ['air_temperature_2m']

        truth = db.get_truth(start_date=20180101, end_date=20180103)
        np.testing.assert_array_equal(truth.indices, np.array([[0, 0], [0, 1], [0, 2], [0, 3], [0, 4], [0, 5], [0, 6], [0, 7], [0, 8]]))


if __name__ == '__main__':
    unittest.main()
