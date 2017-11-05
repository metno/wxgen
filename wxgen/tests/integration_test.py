import os
import unittest
import numpy as np
import tempfile
import netCDF4
import wxgen.driver
import time
np.seterr('raise')


def run_command(command):
   """ Runs a wxgen command line """
   argv = command.split()
   wxgen.driver.main(argv)


def remove(file):
   """ Removes a file """
   os.remove(file)


def file_size(filename):
   """ Returns the number of bytes of a file """
   statinfo = os.stat(filename)
   return statinfo.st_size


def is_valid_file(filename, min_size=3000):
   """ Checks if a file is larger in size than min_size bytes """
   return file_size(filename) > min_size


class IntegrationTest(unittest.TestCase):
   def run_with_image(self, command):
      """
      Runs the wxgen command and appends -f <somefile>.png so that it will write output
      to a temporary png file. Removes the file afterwards.
      """
      fd, imageFile = tempfile.mkstemp(suffix=".png")
      command = command + " -o " + imageFile
      run_command(command)
      self.assertTrue(is_valid_file(imageFile), 3000)
      os.close(fd)
      remove(imageFile)

   def run_with_output(self, command):
      """
      Runs the wxgen command and appends -o <somefile>.nc so that it will write output
      to a temporary nc file.
      """
      fd, file = tempfile.mkstemp(suffix=".nc")
      command = command + " -o " + file
      run_command(command)
      self.assertTrue(is_valid_file(file), 100)
      os.close(fd)
      return file

   def run_with_text(self, command):
      """
      Runs the wxgen command and appends -f <somefile>.txt so that it will write output
      to a temporary txt file. Removes the file afterwards.
      """
      fd, textFile = tempfile.mkstemp(suffix=".txt")
      command = command + " -f " + textFile
      run_command(command)
      self.assertTrue(is_valid_file(textFile, 10))
      os.close(fd)
      remove(textFile)

   def test_valid(self):
      run_command("wxgen")
      run_command("wxgen --version")
      run_command("wxgen sim")
      run_command("wxgen truth")
      run_command("wxgen verif")


class SimTest(IntegrationTest):
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
      self.assertEqual(729, file.dimensions["time"].size)
      self.assertEqual(1, file.dimensions["ensemble_member"].size)
      file.close()

      self.run_with_image("wxgen verif %s %s -m timeseries" % (sim_filename, truth_filename))
      self.run_with_image("wxgen verif %s %s -m variance" % (sim_filename, truth_filename))

      time.sleep(1)

      for filename in [sim_filename, truth_filename]:
         remove(filename)

   def test_stagger(self):
      sim_filename = self.run_with_output("wxgen sim -db examples/database.nc -n 2 -t 730 -g")


class TruthTest(IntegrationTest):
   def test_basic(self):
      self.run_with_output("wxgen truth -db examples/database.nc")
      self.run_with_output("wxgen truth -db examples/database.nc -n 1 -t 365")
      self.run_with_output("wxgen truth -db examples/database.nc -n 2 -t 10")


class VerifTest(IntegrationTest):
   """
   These tests run wxgen on the command-line, but do not test the validity of the output graphics,
   only that they do or do not create errors and that the output image is at least 3 KB in size.
   """

   def test_transform(self):
      sim_filename = self.run_with_output("wxgen sim -db examples/database.nc -n 2 -t 730")
      self.run_with_image("wxgen verif %s -m histogram -tr dryday" % sim_filename)

   def test_aggregator(self):
      sim_filename = self.run_with_output("wxgen sim -db examples/database.nc -n 2 -t 730")
      self.run_with_image("wxgen verif %s -m histogram -a mean" % sim_filename)
      self.run_with_image("wxgen verif %s -m variance -a mean" % sim_filename)
      self.run_with_image("wxgen verif %s -m variance -ea mean" % sim_filename)
      self.run_with_image("wxgen verif %s -m distribution -a mean" % sim_filename)
      self.run_with_image("wxgen verif %s -m distribution -ta median" % sim_filename)

   def test_distribution(self):
      sim_filename = self.run_with_output("wxgen sim -db examples/database.nc -n 2 -t 730")
      self.run_with_image("wxgen verif %s -m distribution" % sim_filename)
      self.run_with_image("wxgen verif %s -m distribution -a mean" % sim_filename)

   def test_timestat(self):
      sim_filename = self.run_with_output("wxgen sim -db examples/database.nc -n 2 -t 730 -v 0")
      self.run_with_image("wxgen verif %s -m timestat" % sim_filename)
      self.run_with_image("wxgen verif %s -m timestat -a variance" % sim_filename)
      self.run_with_image("wxgen verif %s -m timestat -a mean -tr summerday" % sim_filename)
      self.run_with_image("wxgen verif %s -m timestat -a mean -tr summerday -ts 31" % sim_filename)
      self.run_with_image("wxgen verif %s -m timestat -a mean -tm 9" % sim_filename)
      self.run_with_image("wxgen verif %s -m timestat -ea mean -tm 9" % sim_filename)

   def test_legend(self):
      sim_filename1 = self.run_with_output("wxgen sim -db examples/database.nc -n 2 -t 365 -v 0")
      sim_filename2 = self.run_with_output("wxgen sim -db examples/database.nc -n 2 -t 365 -v 0")
      for output in ["variance", "histogram", "timeseries", "timestat"]:
         self.run_with_image("wxgen verif %s %s -m %s -leg 1,2" % (sim_filename1, sim_filename2, output))

   def test_colors_and_styles(self):
      sim_filename1 = self.run_with_output("wxgen sim -db examples/database.nc -n 2 -t 365 -v 0")
      sim_filename2 = self.run_with_output("wxgen sim -db examples/database.nc -n 2 -t 365 -v 0")
      for output in ["variance", "histogram", "timeseries", "timestat"]:
         self.run_with_image("wxgen verif %s %s -m %s -lc red,[1,0.7,0.3] -mfc w,y,k -ls=-,--,- -marker s,None,None" % (sim_filename1, sim_filename2, output))


if __name__ == '__main__':
   unittest.main()
