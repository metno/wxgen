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

if __name__ == '__main__':
   unittest.main()
