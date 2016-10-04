import sys
import argparse
import numpy as np
import matplotlib.pylab as mpl
import wxgen.database
import wxgen.trajectory
import wxgen.generator
import wxgen.verif
import wxgen.output

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
   parser.add_argument('--seed', type=int, default=None, help="Random number seed")

   args = parser.parse_args()

   # Set up database
   if args.seed is not None:
      np.random.seed(args.seed)
   if args.db is None:
      db = wxgen.database.Random(args.n, args.t, args.v)
   else:
      db = wxgen.database.Netcdf(args.db, V=args.v)

   # Generate trajectories
   generator = wxgen.generator.Generator(db)
   initial_state = np.array([275, 0, 0, 0])
   initial_state = np.array([275])
   trajectories = generator.get(args.n, args.t, initial_state)

   # Create output
   if args.type == "plot":
      output = wxgen.output.Timeseries(args.output_filename)
   else:
      output = wxgen.output.Text(args.output_filename)
   output.plot(trajectories)


if __name__ == '__main__':
       main()
