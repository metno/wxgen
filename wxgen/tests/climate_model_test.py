import unittest
import wxgen.climate_model
import numpy as np


class BinTest(unittest.TestCase):
   def test_simple(self):
      model = wxgen.climate_model.Bin(10)
      unixtimes = np.array([wxgen.util.date_to_unixtime(20150101),
         wxgen.util.date_to_unixtime(20150110),
         wxgen.util.date_to_unixtime(20150111)])
      states = model.get(unixtimes)

      self.assertTrue(np.array_equal([[0], [0], [1]], states))

   def test_single(self):
      model = wxgen.climate_model.Bin(10)
      unixtimes = np.array([wxgen.util.date_to_unixtime(20150111)])
      states = model.get(unixtimes)

      self.assertTrue(np.array_equal([[1]], states))


class ComboTest(unittest.TestCase):
   def test_simple(self):
      models = [wxgen.climate_model.Bin(10), wxgen.climate_model.Bin(12)]
      model = wxgen.climate_model.Combo(models)
      unixtimes = np.array([wxgen.util.date_to_unixtime(20150101),
         wxgen.util.date_to_unixtime(20150110),
         wxgen.util.date_to_unixtime(20150111),
         wxgen.util.date_to_unixtime(20150112),
         wxgen.util.date_to_unixtime(20150113)])
      states = model.get(unixtimes)
      self.assertEqual(5, states.shape[0])
      self.assertEqual(2, states.shape[1])
      self.assertTrue(np.array_equal([[0, 0], [0, 0], [1, 0], [1, 0], [1, 1]], states))


if __name__ == '__main__':
   unittest.main()
