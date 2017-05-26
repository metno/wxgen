import unittest
import wxgen.util
import numpy as np


class VariableTest(unittest.TestCase):
   def test_equality(self):
      var1 = wxgen.variable.Variable("Temperature")
      var2 = wxgen.variable.Variable("Temperature")
      self.assertEqual(var1, var2)

   def test_inequality(self):
      var1 = wxgen.variable.Variable("Temperature")
      var2 = wxgen.variable.Variable("Precipitation")
      self.assertTrue(var1 != var2)
      self.assertFalse(var1 == var2)

   def test_pretty(self):
      var1 = wxgen.variable.Variable("Temperature")
      self.assertTrue(var1.pretty() != "")


if __name__ == '__main__':
   unittest.main()
