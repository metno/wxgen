import unittest
import wxgen.util
import numpy as np


class MyTest(unittest.TestCase):
   def test_resize(self):
      vec = np.array([1, 2, 3])
      new = wxgen.util.resize(vec, (3, 4))
      self.assertEqual(new.shape, (3, 4))
      for i in range(0, 4):
         self.assertEqual(new[1, i], 2)

   def test_resize2(self):
      vec = np.array([2.1])
      new = wxgen.util.resize(vec, (3, 4))
      self.assertEqual(new.shape, (3, 4))
      for i in range(0, 4):
         for j in range(0, 3):
            self.assertEqual(new[j, i], 2.1)

   def test_resize3(self):
      vec = np.array([[9, 6, 2], [4, 1, 3]])  # size: (2,3)
      new = wxgen.util.resize(vec, (4, 9))
      self.assertEqual(new.shape, (4, 9))
      self.assertTrue(np.array_equal(new[:, 0], np.array([9, 4, 9, 4])))
      self.assertTrue(np.array_equal(new[:, 4], np.array([6, 1, 6, 1])))


class TestRandomWeighted(unittest.TestCase):
   def test_simple(self):
      weights = np.array([1, 2, 3])
      for i in range(0, 5):
         I = wxgen.util.random_weighted(weights, "top2")
         self.assertTrue(I >= 1 and I <= 2)
      weights = np.array([3, 1, 2])
      I = wxgen.util.random_weighted(weights, "top1")
      self.assertTrue(I == 0)

   def test_top5_ties(self):
      weights = np.array([0, 3, 1, 3, 1, 3, 1, 3, 1, 3])
      I = wxgen.util.random_weighted(weights, "top2")


class TestClimatology(unittest.TestCase):
   array = np.array([[1, 2], [2, 2], [3, 1], [12, 0], [-2, 2]])  # 5 times, 2 members

   def test_simple(self):
      clim = wxgen.util.climatology(self.array)
      self.assertEqual(1, len(clim.shape))
      self.assertEqual(5, len(clim))
      self.assertTrue(np.array_equal(clim, np.array([1.5, 2, 2, 6, 0])))

   def test_window(self):
      clim = wxgen.util.climatology(self.array, window=3)
      self.assertEqual(1, len(clim.shape))
      self.assertEqual(5, len(clim))
      self.assertTrue(np.min(np.isclose(clim, np.array([5.0/3, 11.0/6, 10.0/3, 8.0/3, 2]))))


if __name__ == '__main__':
   unittest.main()
