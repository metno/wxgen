import sys
import argparse
import numpy as np
import matplotlib.pylab as mpl
import wxgen.database
import wxgen.trajectory
import wxgen.generator
import wxgen.metric
import wxgen.verif
import wxgen.output
import wxgen.version

#@profile
def run(argv):
   if 0:
      if len(sys.argv) < 3:
         print "Weather generator"
         print "usage: wxgen N T [db]"
         print ""
         print "Arguments:"
         print "  N: Number of trajectories"
         print "  T: Number of days in trajectory"
         print "  db: Filename of Netcdf database of trajectory"

         sys.exit()

   parser = argparse.ArgumentParser(description="Weather generator")
   parser.add_argument('-n', type=int, help="Number of trajectories", required=True)
   parser.add_argument('-t', type=int, help="Length of trajectory", required=True)
   parser.add_argument('-v', type=int, help="Number of variables", required=False)
   parser.add_argument('--type', type=str, default="plot", help="Output type (text or plot)")
   parser.add_argument('--db', type=str, default=None, help="Filename of NetCDF database")
   parser.add_argument('-o', type=str, default=None, help="Output filename", dest="output_filename")
   parser.add_argument('-m', type=str, default="rmsd", help="Metric for matching states (currently only rmsd)")
   parser.add_argument('--seed', type=int, default=None, help="Random number seed")
   parser.add_argument('--debug', help="Display debug information", action="store_true")
   parser.add_argument('--version', action="version", version=wxgen.version.__version__)

   args = parser.parse_args()

   # Set up database
   if args.seed is not None:
      np.random.seed(args.seed)
   if args.db is None:
      # Don't use args.t as the segment length, because then you never get to join
      # Don't use args.n as the number of segments, because then you never get to join
      db = wxgen.database.Random(100, 10, args.v)
   else:
      db = wxgen.database.Netcdf(args.db, V=args.v)
   if args.debug:
      db.info()

   # Generate trajectories
   metric = wxgen.metric.get(args.m)
   generator = wxgen.generator.Generator(db, metric)
   initial_state = np.array([270, 0, 0, 0])
   initial_state = np.array([270])
   trajectories = generator.get(args.n, args.t, initial_state)

   # Create output
   if args.type == "plot":
      output = wxgen.output.Timeseries(db, args.output_filename)
   elif args.type == "db":
      output = wxgen.output.Database(db, args.output_filename)
   else:
      output = wxgen.output.Text(db, args.output_filename)
   output.plot(trajectories)


if __name__ == '__main__':
       main()
