import unittest
import wxgen.generator
import wxgen.database
import numpy as np


class MyTest(unittest.TestCase):
   def test_getOffsets(self):
      V = 3
      N = 4
      T = 20
      db = wxgen.database.Random(N, T, V)
      generator = wxgen.generator.Generator(db)
      traj = generator.get(N, T)
      self.assertEqual(N, len(traj))
      self.assertEqual(T, traj[0].shape[0])
      self.assertEqual(V, traj[0].shape[1])


if __name__ == '__main__':
   unittest.main()
