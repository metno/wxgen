import os
import unittest
import numpy as np
import tempfile
import netCDF4
import wxgen.driver
np.seterr('raise')


class IntegrationTest(unittest.TestCase):
   """
   These tests run wxgen on the command-line, but do not test the validity of the
   output graphics, only that they do or do not create errors.
   """

   @staticmethod
   def run_command(command):
      """ Runs a wxgen command line """
      argv = command.split()
      wxgen.driver.main(argv)

   @staticmethod
   def remove(file):
      """ Removes a file """
      os.remove(file)

   @staticmethod
   def file_size(filename):
      """ Returns the number of bytes of a file """
      statinfo = os.stat(filename)
      return statinfo.st_size

   @staticmethod
   def is_valid_file(filename, min_size=3000):
      """ Checks if a file is larger in size than min_size bytes """
      return IntegrationTest.file_size(filename) > min_size

   def run_with_image(self, command):
      """
      Runs the wxgen command and appends -f <somefile>.png so that it will write output
      to a temporary png file. Removes the file afterwards.
      """
      fd, imageFile = tempfile.mkstemp(suffix=".png")
      command = command + " -o " + imageFile
      self.run_command(command)
      self.assertTrue(self.is_valid_file(imageFile), 3000)
      os.close(fd)
      self.remove(imageFile)

   def run_with_output(self, command):
      """
      Runs the wxgen command and appends -o <somefile>.nc so that it will write output
      to a temporary nc file.
      """
      fd, file = tempfile.mkstemp(suffix=".nc")
      command = command + " -o " + file
      self.run_command(command)
      os.close(fd)
      return file

   def run_with_text(self, command):
      """
      Runs the wxgen command and appends -f <somefile>.txt so that it will write output
      to a temporary txt file. Removes the file afterwards.
      """
      fd, textFile = tempfile.mkstemp(suffix=".txt")
      command = command + " -f " + textFile
      self.run_command(command)
      self.assertTrue(self.is_valid_file(textFile, 10))
      os.close(fd)
      self.remove(textFile)

   def valid(self):
      self.run_command("wxgen")
      self.run_command("wxgen --version")
      self.run_command("wxgen sim")
      self.run_command("wxgen truth")
      self.run_command("wxgen verif")

   def test_README(self):
      sim_filename = self.run_with_output("wxgen sim -db examples/database.nc -n 10 -t 100")
      file = netCDF4.Dataset(sim_filename, 'r')
      self.assertTrue("time" in file.dimensions)
      self.assertTrue("ensemble_member" in file.dimensions)
      self.assertEqual(100, file.dimensions["time"].size)
      self.assertEqual(10, file.dimensions["ensemble_member"].size)
      file.close()

      truth_filename = self.run_with_output("wxgen truth -db examples/database.nc")
      file = netCDF4.Dataset(truth_filename, 'r')
      self.assertTrue("time" in file.dimensions)
      self.assertTrue("ensemble_member" in file.dimensions)
      self.assertEqual(364, file.dimensions["time"].size)
      self.assertEqual(1, file.dimensions["ensemble_member"].size)
      file.close()

      self.run_with_image("wxgen verif %s -truth %s -m timeseries" % (sim_filename, truth_filename))
      self.run_with_image("wxgen verif %s -truth %s -m variance" % (sim_filename, truth_filename))

      for filename in [sim_filename, truth_filename]:
         self.remove(filename)


if __name__ == '__main__':
   unittest.main()
