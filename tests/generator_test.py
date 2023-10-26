import unittest
import wxgen.generator
import wxgen.database
import numpy as np


class GeneratorTest(unittest.TestCase):
    def test_get(self):
        V = 3
        N = 4
        T = 20
        db = wxgen.database.Random(N, T, V)
        generator = wxgen.generator.Generator(db)
        traj = generator.get(N, T)
        self.assertEqual(N, len(traj))
        self.assertEqual(T, traj[0].length)
        self.assertEqual(2, traj[0].indices.shape[1])


if __name__ == '__main__':
    unittest.main()
