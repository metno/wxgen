import sys
import argparse
import numpy as np
import wxgen.climate_model
import wxgen.database
import wxgen.trajectory
import wxgen.generator
import wxgen.metric
import wxgen.output
import wxgen.plot
import wxgen.version


def main(argv):
   parser = argparse.ArgumentParser(prog="wxgen", description="Hybrid weather generator, combining stochastic and physical modelling")
   parser.add_argument('--version', action="version", version=wxgen.version.__version__)
   subparsers = parser.add_subparsers(title="Choose one of these commands", dest="command")

   """
   Simulation driver
   """
   sp = dict()
   sp["sim"] = subparsers.add_parser('sim', help='Create simulated scenarios')
   sp["sim"].add_argument('-n', metavar="NUM", type=int, help="Number of trajectories", required=True)
   sp["sim"].add_argument('-t', metavar="DAYS", type=int, help="Length of trajectory", required=True)
   sp["sim"].add_argument('-m', default="rmsd", help="Metric for matching states (currently only rmsd)", dest="metric")
   sp["sim"].add_argument('-rs', type=int, help="Random number seed", dest="seed")
   sp["sim"].add_argument('-w', help="Weights for each variable when joining (comma-separated)", dest="weights")
   sp["sim"].add_argument('-i', help="Initial state", dest="initial")
   sp["sim"].add_argument('-j', type=int, metavar="NUM", help="How many times should segments be prejoined?", dest="prejoin")
   sp["sim"].add_argument('-b', type=int, metavar="DAYS", help="Length of database bins", dest="bin_width")

   """
   Truth trajetory driver
   """
   sp["truth"] = subparsers.add_parser('truth', help='Create truth scenario')
   sp["truth"].add_argument('-ed', metavar="YYYYMMDD", type=int, help="End date of trajectory", dest="end_date")
   sp["truth"].add_argument('-sd', metavar="YYYYMMDD", type=int, help="Start date of trajectory", dest="start_date")

   for driver in ["sim", "truth"]:
      sp[driver].add_argument('-s', default="agg", help="Output scale (agg, large, small)", dest="scale")
      sp[driver].add_argument('-db', metavar="FILENAME", help="Filename of NetCDF database")
      sp[driver].add_argument('-dbtype', metavar="TYPE", help="Database type (netcdf, random, lorenz63). If --db is provided, then --db_type is automatically set to 'netcdf'. If neither --db nor --db_type is set, then --db_type is automatically set to 'random'.")
      sp[driver].add_argument('-o', metavar="FILENAME", help="Output filename", dest="filename", required=True)

   """
   Verification driver
   """
   sp["verif"] = subparsers.add_parser('verif', help='Verify trajectories')
   sp["verif"].add_argument('files', help="Input files", nargs="*")
   sp["verif"].add_argument('-fs', type=wxgen.util.parse_ints, default=[10, 5], help="Figure size: width,height")
   sp["verif"].add_argument('-m', metavar="METRIC", help="Verification metric", dest="metric")
   sp["verif"].add_argument('-o', metavar="FILENAME", help="Output filename", dest="filename")
   sp["verif"].add_argument('-truth', metavar="FILENAME", help="File with truth scenario", dest="truth")
   sp["verif"].add_argument('-xlim', type=wxgen.util.parse_ints, help="x-axis limits: lower,upper")
   sp["verif"].add_argument('-xlog', help="X-axis log scale", action="store_true")
   sp["verif"].add_argument('-ylim', help="y-axis limits: lower,upper")
   sp["verif"].add_argument('-ylog', help="Y-axis log scale", action="store_true")
   sp["verif"].add_argument('-r', dest="thresholds", help="Thresholds for use in plots", required=False, type=wxgen.util.parse_numbers)

   """
   Common options
   """
   for driver in sp.keys():
      sp[driver].add_argument('-v', metavar="INDICES", help="Which variables to use? Use indices, starting at 0.", required=False, type=wxgen.util.parse_ints, dest="vars")
      sp[driver].add_argument('--debug', help="Display debug information", action="store_true")

   if len(argv) < 2:
      parser.print_help()
      return
   elif len(argv) == 2 and argv[1] in sp.keys():
      sp[argv[1]].print_help()
      return

   args = parser.parse_args(argv[1:])

   """
   Run commands
   """
   wxgen.util.DEBUG = args.debug
   if args.command == "sim":
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
      model = get_climate_model(args)
      generator = wxgen.generator.LargeScale(db, metric, model=model)
      generator.prejoin = args.prejoin
      trajectories = generator.get(args.n, args.t, initial_state)

      # Create output
      output = wxgen.output.Netcdf(args.filename)
      output.write(trajectories, db, args.scale)

   elif args.command == "truth":
      db = get_db(args)
      trajectory = db.get_truth(args.start_date, args.end_date)
      output = wxgen.output.Netcdf(args.filename)
      output.write([trajectory], db, args.scale)

   elif args.command == "verif":
      plot = wxgen.plot.get(args.metric)()
      plot.filename = args.filename
      plot.xlim = args.xlim
      plot.xlog = args.xlog
      plot.ylim = args.ylim
      plot.ylog = args.ylog
      plot.vars = args.vars
      plot.thresholds = args.thresholds
      plot.fig_size = args.fs
      truth = None
      sims = None
      if args.truth is not None:
         truth = wxgen.database.Netcdf(args.truth, None)
      if args.files is not None:
         sims = [wxgen.database.Netcdf(file, None) for file in args.files]
      plot.plot(sims, truth)


def get_db(args):
   # Set up database
   if args.command == "sim" and args.seed is not None:
      np.random.seed(args.seed)

   model = get_climate_model(args)

   if args.db is None:
      # Don't use args.t as the segment length, because then you never get to join
      # Don't use args.n as the number of segments, because then you never get to join
      if args.db_type is None or args.db_type == "random":
         db = wxgen.database.Random(100, 10, 3, model=model)
      elif args.db_type == "lorenz63":
         db = wxgen.database.Lorenz63(10, 500, model=model)
      else:
         wxgen.util.error("Cannot understand --db_type %s" % args.db_type)
   else:
      db = wxgen.database.Netcdf(args.db, args.vars, model=model)

   if args.debug:
      db.info()

   return db


def get_metric(args):
   if args.weights is not None:
      weights = np.array(wxgen.util.parse_numbers(args.weights))
      # if len(weights) != len(db.variables):
      #    wxgen.util.error("Weights must match the number of variables (%d)" % (V))

      if args.metric == "rmsd":
         metric = wxgen.metric.Rmsd(weights)
      elif args.metric == "exp":
         metric = wxgen.metric.Exp(weights)
   else:
      metric = wxgen.metric.get(args.metric)

   return metric


def get_climate_model(args):
   model = None
   try:
      # args might not have bin_width
      if args.bin_width is not None:
         model = wxgen.climate_model.Bin(args.bin_width)
   except:
      pass
   return model


if __name__ == '__main__':
   main()
