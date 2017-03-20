import unittest
import wxgen.util
import wxgen.aggregator
import numpy as np


class AggregatorTest(unittest.TestCase):
   def test_consecutive(self):
      agg = wxgen.aggregator.Consecutive()
      self.assertEqual(agg(np.array([0, 0, 1, 1, 1, 0, 0, 0])), 3)
      self.assertEqual(agg(np.array([0, 0, 0, 0, 1, 0, 0, 0])), 4)
      self.assertEqual(agg(np.array([1, 1, 1])), 0)
      self.assertEqual(agg(np.array([0, 0])), 2)
      self.assertEqual(agg(np.array([1, 0, 1, 0, 1, 0, 0])), 2)
      self.assertEqual(agg(np.array([1, 0, 1, 0, 1, 0])), 1)
      self.assertEqual(agg(np.array([0])), 1)
      self.assertEqual(agg(np.array([1])), 0)


if __name__ == '__main__':
   unittest.main()
