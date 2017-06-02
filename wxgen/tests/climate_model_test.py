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

   def test_last_bin(self):
      def get_state(date):
         unixtimes = np.array([wxgen.util.date_to_unixtime(date)])
         states = model.get(unixtimes)
         return states[0,0]

      # Test that we do not get a small bin at the end of the year
      # That is, 20151230 should be assigned to bin 11 not bin 12
      model = wxgen.climate_model.Bin(30)
      self.assertEqual(11, get_state(20151210))
      self.assertEqual(11, get_state(20151220))
      self.assertEqual(11, get_state(20151230))

      # The last bin has less than half of the first bin, therefore
      # don't create the last bin.
      model = wxgen.climate_model.Bin(250)
      self.assertEqual(0, get_state(20150101))
      self.assertEqual(0, get_state(20150601))
      self.assertEqual(0, get_state(20151201))

      # The last bin has more than half of the first bin, therefore
      # do create the last bin.
      model = wxgen.climate_model.Bin(200)
      self.assertEqual(0, get_state(20150101))
      self.assertEqual(0, get_state(20150601))
      self.assertEqual(1, get_state(20151201))

      model = wxgen.climate_model.Bin(1000)
      self.assertEqual(0, get_state(20150101))
      self.assertEqual(0, get_state(20150601))
      self.assertEqual(0, get_state(20151201))


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
