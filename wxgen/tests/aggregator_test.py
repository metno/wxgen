import unittest
import wxgen.util
import wxgen.aggregator
import numpy as np


class AggregatorTest(unittest.TestCase):
   def test_consecutive(self):
      agg = wxgen.aggregator.Consecutive()
      self.assertEqual(agg(np.array([1, 1, 0, 0, 0, 1, 1, 1])), 3)
      self.assertEqual(agg(np.array([1, 1, 1, 1, 0, 1, 1, 1])), 4)
      self.assertEqual(agg(np.array([0, 0, 0])), 0)
      self.assertEqual(agg(np.array([1, 1])), 2)
      self.assertEqual(agg(np.array([0, 1, 0, 1, 0, 1, 1])), 2)
      self.assertEqual(agg(np.array([0, 1, 0, 1, 0, 1])), 1)
      self.assertEqual(agg(np.array([1])), 1)
      self.assertEqual(agg(np.array([0])), 0)

   def test_consecutive_axis(self):
      agg = wxgen.aggregator.Consecutive()
      # self.assertEqual(agg(np.array([1, 1, 0, 0, 0, 1, 1, 1]), axis=0), 3)
      ar = np.array([[1, 0, 0], [1, 1, 0]])
      self.assertTrue(np.array_equal(agg(ar, axis=0), [2, 1, 0]))
      self.assertTrue(np.array_equal(agg(ar, axis=1), [1, 2]))
      np.array_equal(agg(np.array([1, 1, 0]), axis=0), [1, 1, 0])


if __name__ == '__main__':
   unittest.main()
