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
   parser = argparse.ArgumentParser(description="Weather generator")
   parser.add_argument('-n', type=int, help="Number of trajectories", required=True)
   parser.add_argument('-t', type=int, help="Length of trajectory", required=True)
   parser.add_argument('-v', type=int, help="Number of variables", required=False)
   parser.add_argument('--which_vars', type=str, help="Which variables to plot?", required=False)
   parser.add_argument('--type', type=str, default="timeseries", help="Output type (text, netcdf, or plot)")
   parser.add_argument('--db', type=str, default=None, help="Filename of NetCDF database")
   parser.add_argument('--db_type', type=str, default=None, help="Database type (netcdf, random, lorenz63). If --db is provided, then --db_type is automatically set to 'netcdf'. If neither --db nor --db_type is set, then --db_type is automatically set to 'random'.")
   parser.add_argument('-o', type=str, default=None, help="Output filename", dest="output_filename")
   parser.add_argument('-m', type=str, default="rmsd", help="Metric for matching states (currently only rmsd)")
   parser.add_argument('-fs', type=str, default=None, help="Figure size: width,height")
   parser.add_argument('--xlog', help="x-axis limits: lower,upper", action="store_true")
   parser.add_argument('--ylog', help="y-axis limits: lower,upper", action="store_true")
   parser.add_argument('--xlim', type=str, default=None, help="x-axis limits: lower,upper")
   parser.add_argument('--ylim', type=str, default=None, help="y-axis limits: lower,upper")
   parser.add_argument('--seed', type=int, default=None, help="Random number seed")
   parser.add_argument('--debug', help="Display debug information", action="store_true")
   parser.add_argument('--version', action="version", version=wxgen.version.__version__)
   parser.add_argument('--weights', type=str)
   parser.add_argument('--initial', type=str, default=None, help="Initial state")

   args = parser.parse_args()

   # Set up database
   if args.seed is not None:
      np.random.seed(args.seed)
   if args.db is None:
      # Don't use args.t as the segment length, because then you never get to join
      # Don't use args.n as the number of segments, because then you never get to join
      if args.db_type is None or args.db_type == "random":
         db = wxgen.database.Random(100, 10, args.v)
      elif args.db_type == "lorenz63":
         db = wxgen.database.Lorenz63(10, 500)
      else:
         wxgen.util.error("Cannot understand --db_type %s" % args.db_type)
   else:
      db = wxgen.database.Netcdf(args.db, V=args.v)
   if args.debug:
      db.info()
   V = len(db.variables)

   if args.initial is None:
      initial_state = None
   else:
      initial_state = np.array(wxgen.util.parse_numbers(args.initial))
      if len(initial_state) != V:
         wxgen.util.error("Initial state must match the number of variables (%d)" % (V))

   # Generate trajectories
   if args.weights is not None:
      weights = np.array(wxgen.util.parse_numbers(args.weights))
      if len(weights) != len(db.variables):
         wxgen.util.error("Weights must match the number of variables (%d)" % (V))

      if args.m == "rmsd":
         metric = wxgen.metric.Rmsd(weights)
      elif args.m == "exp":
         metric = wxgen.metric.Exp(weights)
   else:
      metric = wxgen.metric.get(args.m)

   generator = wxgen.generator.Generator(db, metric)
   trajectories = generator.get(args.n, args.t, initial_state)

   # Create output
   output = wxgen.output.get(args.type)(db)
   output.filename = args.output_filename
   output.fig_size = wxgen.util.parse_numbers(args.fs)
   output.xlim = wxgen.util.parse_numbers(args.xlim)
   output.ylim = wxgen.util.parse_numbers(args.ylim)
   output.xlog = args.xlog
   output.ylog = args.ylog
   output.which_vars = wxgen.util.parse_numbers(args.which_vars)
   print output.which_vars
   output.plot(trajectories)


if __name__ == '__main__':
       main()
