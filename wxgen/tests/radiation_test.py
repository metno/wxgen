import unittest
import wxgen.radiation
import numpy as np


class SwingTest(unittest.TestCase):
   def test(self):
      for i in range(24):
         jday = np.array([90, 90])
         hour = np.array([i, i])
         lat = np.array([60, 60])
         lon = np.array([10, 10])
         cloud_cover = np.array([0.1, 0.1])
         pressure = np.array([1013.25, 1013.25])
         temperature = np.array([0, 30])
         sunrad = wxgen.radiation.swing(jday, hour, lat, lon, cloud_cover, pressure, temperature)
         # self.assertEqual(312, sunrad[0])


if __name__ == '__main__':
   unittest.main()
