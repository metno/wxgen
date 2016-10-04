import sys
import matplotlib.pylab as mpl
import wxgen.database
import wxgen.trajectory
import wxgen.generator
import wxgen.verif
import wxgen.output
import numpy as np

#@profile
def run(argv):
   if len(sys.argv) < 3:
      print "Weather generator"
      print "usage: wxgen N T [db]"
      print ""
      print "Arguments:"
      print "  N: Number of trajectories"
      print "  T: Number of days in trajectory"
      print "  db: Filename of Netcdf database of trajectory"

      sys.exit()

   N = int(sys.argv[1])
   T = int(sys.argv[2])
   np.random.seed(1)
   if len(sys.argv) == 3:
      V = 2
      db = wxgen.database.Random(N, T, V)
   else:
      db_filename = sys.argv[3]
      db = wxgen.database.Netcdf(db_filename)

   generator = wxgen.generator.Generator(db)
   initial_state = np.array([275, 0, 0, 0])
   initial_state = np.array([275])
   trajectories = generator.get(N, T, initial_state)

   #output = wxgen.output.Text("test.txt")
   output = wxgen.output.Timeseries()
   output.plot(trajectories)


if __name__ == '__main__':
       main()
