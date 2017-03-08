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
import wxgen.plot
import wxgen.version

def main(argv):
   parser = argparse.ArgumentParser(description="Hybrid weather generator, combining stochastic and physical modelling")
   parser.add_argument('--version', action="version", version=wxgen.version.__version__)
   subparsers = parser.add_subparsers(title="Choose one of these commands", dest="command")

   """
   Simulation driver
   """
   sp = dict()
   sp["sim"] = subparsers.add_parser('sim', help='Create simulated scenarios')
   sp["sim"].add_argument('-n', type=int, help="Number of trajectories", required=True)
   sp["sim"].add_argument('-t', type=int, help="Length of trajectory", required=True)
   sp["sim"].add_argument('--vars', metavar="INDICES", help="Which variables to use? Use indices, starting at 0.", required=False, type=wxgen.util.parse_ints)
   sp["sim"].add_argument('--db', type=str, default=None, help="Filename of NetCDF database")
   sp["sim"].add_argument('--db_type', type=str, default=None, help="Database type (netcdf, random, lorenz63). If --db is provided, then --db_type is automatically set to 'netcdf'. If neither --db nor --db_type is set, then --db_type is automatically set to 'random'.")
   sp["sim"].add_argument('-o', metavar="FILENAME", help="Output filename", dest="filename")
   sp["sim"].add_argument('--scale', type=str, default="agg", help="Output scale (agg, large, small)")
   sp["sim"].add_argument('-m', type=str, default="rmsd", help="Metric for matching states (currently only rmsd)")
   sp["sim"].add_argument('--seed', type=int, default=None, help="Random number seed")
   sp["sim"].add_argument('--debug', help="Display debug information", action="store_true")
   sp["sim"].add_argument('--weights', type=str)
   sp["sim"].add_argument('--initial', type=str, default=None, help="Initial state")

   # p_sim.add_argument('--type', type=str, default="timeseries", help="Output type (text, netcdf, or plot)")
   # p_sim.add_argument('-fs', type=str, default=None, help="Figure size: width,height")
   # p_sim.add_argument('--xlog', help="x-axis limits: lower,upper", action="store_true")
   # p_sim.add_argument('--ylog', help="y-axis limits: lower,upper", action="store_true")
   # p_sim.add_argument('--xlim', type=str, default=None, help="x-axis limits: lower,upper")
   # p_sim.add_argument('--ylim', type=str, default=None, help="y-axis limits: lower,upper")

   """
   Truth trajetory driver
   """
   sp["truth"] = subparsers.add_parser('truth', help='Create truth scenario')
   sp["truth"].add_argument('-o', metavar="FILENAME", help="Output filename", dest="filename", required=True)
   sp["truth"].add_argument('-d', type=wxgen.util.parse_dates, help="Limit trajectory to these dates (e.g. 20010101:20031231)", dest="dates")
   sp["truth"].add_argument('--db', metavar="FILENAME", help="Filename of NetCDF database")
   sp["truth"].add_argument('--db_type', metavar="TYPE", help="Database type (netcdf, random, lorenz63). If --db is provided, then --db_type is automatically set to 'netcdf'. If neither --db nor --db_type is set, then --db_type is automatically set to 'random'.")
   sp["truth"].add_argument('--vars', metavar="INDICES", help="Which variables to use? Use indices, starting at 0.", required=False, type=wxgen.util.parse_ints)
   sp["truth"].add_argument('--debug', help="Display debug information", action="store_true")

   """
   Verification driver
   """
   sp["verif"] = subparsers.add_parser('verif', help='Verify trajectories')
   sp["verif"].add_argument('-n', type=int, help="Number of trajectories")
   sp["verif"].add_argument('files', help="Input files", nargs="+")

   if len(sys.argv) < 2:
      parser.print_help()
      sys.exit(1)
   elif len(sys.argv) == 2 and sys.argv[1] in sp.keys():
      sp[sys.argv[1]].print_help()
      sys.exit(1)

   args = parser.parse_args()


   """
   Run commands
   """
   if args.command == "sim":
      wxgen.util.DEBUG = args.debug
      db = get_db(args)
      V = len(db.variables)
      if args.initial is None:
         initial_state = None
      else:
         initial_state = np.array(wxgen.util.parse_numbers(args.initial))
         if len(initial_state) != V:
            wxgen.util.error("Initial state must match the number of variables (%d)" % (V))

      # Generate trajectories
      metric = get_metric(args)
      generator = wxgen.generator.LargeScale(db, metric)
      trajectories = generator.get(args.n, args.t, initial_state)

      # Create output
      output = wxgen.output.Netcdf(args.filename)
      output.write(trajectories, db, args.scale)

   elif args.command == "truth":
      db = get_db(args)
      trajectory = db.get_truth()
      output = wxgen.output.Netcdf(args.filename)
      output.write([trajectory], db, "agg")

   elif args.command == "verif":
      plot = wxgen.plot.Variance()
      db = wxgen.database.Netcdf(args.files[0], None)
      trajectories = [db.get(i) for i in range(db.num)]
      plot.plot(trajectories, db)
      pass


def get_db(args):
   # Set up database
   if args.command == "sim" and args.seed is not None:
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
      db = wxgen.database.Netcdf(args.db, args.vars)

   if args.debug:
      db.info()

   return db


def get_metric(args):
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


   return metric


if __name__ == '__main__':
   main()
