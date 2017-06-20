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
   parser, sp = get_parsers()

   # Show help message when no options are provided
   if len(argv) < 2:
      parser.print_help()
      return
   elif len(argv) == 2 and argv[1] == "--version":
      print(wxgen.version.__version__)
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
      metric = get_metric(args, db)
      generator = wxgen.generator.LargeScale(db, metric)
      generator.prejoin = args.prejoin
      generator.policy = args.policy
      trajectories = generator.get(args.n, args.t, initial_state)

      # Create output
      output = wxgen.output.Netcdf(args.filename)
      output.lat = args.lat
      output.lon = args.lon
      output.write(trajectories, db, args.scale)

   elif args.command == "truth":
      db = get_db(args)
      if args.n is None and args.t is None:
        trajectories = [db.get_truth(args.start_date, args.end_date)]
      else:
         """
         Create a number of ensemble members, by sampling long timeseries from database (no
         joining). This is determined by -sd, -ed, -n, and -t. -sd -ed sets which dates are allowed
         to use from the database. Then -n sets the number of scenarios.

         We only want to create scenarios that all start at the same time of the year.
         """
         if args.n is None:
            wxgen.util.error("-n not specified")
         if args.t is None:
            wxgen.util.error("-t not specified")

         # Determine allowable dates from database
         start_date = args.start_date
         end_date = args.end_date
         if args.start_date is None:
            start_date = wxgen.util.unixtime_to_date(np.min(db.inittimes))
         if args.end_date is None:
            end_date = wxgen.util.unixtime_to_date(np.max(db.inittimes))
         dates = np.array(wxgen.util.parse_dates("%d:%d" % (start_date, end_date)))

         # Figure out which time indices are possible starting dates, by finding dates that have the
         # same day of year as the first date of the allowable dates
         months = dates / 100 % 100
         days = dates % 100
         start_day = start_date % 100
         start_month = start_date / 100 % 100
         Ipossible_start_days = np.where((months == start_month) & (days == start_day))[0]
         if args.n > len(Ipossible_start_days):
            wxgen.util.warning("Not enough possible starting days (%d < %d)" % (len(Ipossible_start_days), args.n))

         trajectories = list()
         for n in range(min(args.n, len(Ipossible_start_days))):
            s = dates[Ipossible_start_days[n]]
            e = wxgen.util.get_date(s, args.t)
            if e < end_date:
               wxgen.util.debug("Member %d dates: %d - %d" % (n, s, e))
               trajectory = db.get_truth(s, e)
               trajectories += [trajectory]
            else:
               wxgen.util.debug("Skipping member %d: Goes outside date range" % n, "yellow")

      output = wxgen.output.Netcdf(args.filename)
      output.lat = args.lat
      output.lon = args.lon
      output.write(trajectories, db, args.scale)

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
      plot.clim = args.clim
      plot.cmap = args.cmap
      plot.lat = args.lat
      plot.lon = args.lon
      plot.timemod = args.timemod
      plot.timescale = args.timescale
      plot.scale = args.scale

      transform = get_transform(args)
      if transform is not None:
         if not plot.supports_transform:
            wxgen.util.error("Plot does not support -tr")
         plot.transform = transform
      aggregator = get_aggregator(args)
      if aggregator is not None:
         if not plot.supports_aggregator:
            wxgen.util.error("Plot does not support -a")
         plot.aggregator = aggregator

      sims = [wxgen.database.Netcdf(file, None) for file in args.files]
      plot.plot(sims)


def get_parsers():
   """ Sets up the argparser object for wxgen

   Returns:
      parser: The main argparse object
      sp: A list of subparsers
   """
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
   sp["sim"].add_argument('-m', default="rmsd", help="Metric for matching states", dest="metric", choices=get_module_names(wxgen.metric))
   sp["sim"].add_argument('-rs', type=int, help="Random number seed", dest="seed")
   sp["sim"].add_argument('-w', help="Weights for each variable when joining (comma-separated)", dest="weights")
   sp["sim"].add_argument('-i', help="Initial state", dest="initial")
   sp["sim"].add_argument('-j', type=int, metavar="NUM", help="How many times should segments be prejoined?", dest="prejoin")
   sp["sim"].add_argument('-b', type=int, metavar="DAYS", help="Length of database bins", dest="bin_width")
   sp["sim"].add_argument('-p', default="top5", metavar="POLICY", help="Randomization policy. One of 'random', 'top<N>'", dest="policy")

   """
   Truth trajetory driver
   """
   sp["truth"] = subparsers.add_parser('truth', help='Create truth scenario')
   sp["truth"].add_argument('-sd', metavar="YYYYMMDD", type=int, help="Earliest date to use from database", dest="start_date")
   sp["truth"].add_argument('-ed', metavar="YYYYMMDD", type=int, help="Latest date to use from database", dest="end_date")
   sp["truth"].add_argument('-n', metavar="NUM", type=int, help="Number of trajectories (if -n and -t are unspecified, create one trajectory with all data)")
   sp["truth"].add_argument('-t', metavar="DAYS", type=int, help="Length of trajectory")

   for driver in ["sim", "truth"]:
      sp[driver].add_argument('-db', metavar="FILENAME", help="Filename of NetCDF database")
      sp[driver].add_argument('-dbtype', metavar="TYPE", help="Database type (netcdf, random, lorenz63). If --db is provided, then --dbtype is automatically set to 'netcdf'. If neither --db nor --dbtype is set, then --dbtype is automatically set to 'random'.")
      sp[driver].add_argument('-o', metavar="FILENAME", help="Output filename", dest="filename", required=True)
      sp[driver].add_argument('-wl', type=int, default=0, metavar="NUM", help="Number of wavelet levels.  If 0 (default), don't use wavelets.", dest="wavelet_levels")

   """
   Verification driver
   """
   sp["verif"] = subparsers.add_parser('verif', help='Verify trajectories')
   sp["verif"].add_argument('files', help="Input files", nargs="+")
   sp["verif"].add_argument('-fs', type=wxgen.util.parse_ints, default=[10, 5], help="Figure size: width,height")
   sp["verif"].add_argument('-m', help="Verification metric", dest="metric", required=True, choices=get_module_names(wxgen.plot))
   sp["verif"].add_argument('-o', metavar="FILENAME", help="Output filename", dest="filename")
   sp["verif"].add_argument('-xlim', type=wxgen.util.parse_ints, help="x-axis limits: lower,upper")
   sp["verif"].add_argument('-xlog', help="X-axis log scale", action="store_true")
   sp["verif"].add_argument('-ylim', type=wxgen.util.parse_ints, help="y-axis limits: lower,upper")
   sp["verif"].add_argument('-ylog', help="Y-axis log scale", action="store_true")
   sp["verif"].add_argument('-r', dest="thresholds", help="Thresholds for use in plots", required=False, type=wxgen.util.parse_numbers)
   sp["verif"].add_argument('-tr', dest="transform", help="Transform for use in plots", choices=get_module_names(wxgen.transform))
   aggregators = [mod for mod in get_module_names(wxgen.aggregator) if mod != "quantile"]
   sp["verif"].add_argument('-a', dest="aggregator", help="Aggregator for use in plots. One of: " + ', '.join(aggregators) + " or a number between 0 and 1 representing a quantile")
   sp["verif"].add_argument('-clim', type=wxgen.util.parse_numbers, help="Colorbar limits (lower,upper)")
   sp["verif"].add_argument('-cmap', help="Colormap (e.g. jet, RdBu, Blues_r)")
   sp["verif"].add_argument('-tm', type=int, help="Time modulus (in days)", dest="timemod")
   sp["verif"].add_argument('-ts', default=1, type=int, help="Time scale (in days)", dest="timescale")

   """
   Common options
   """
   for driver in sp.keys():
      sp[driver].add_argument('-v', metavar="INDICES", help="Which variables to use? Use indices, starting at 0.", required=False, type=wxgen.util.parse_ints, dest="vars")
      sp[driver].add_argument('--debug', help="Display debug information", action="store_true")
      sp[driver].add_argument('-s', default="large", help="Output scale", choices=["agg", "large", "small"], dest="scale")
      sp[driver].add_argument('-lat', type=float, help="Lookup latitude")
      sp[driver].add_argument('-lon', type=float, help="Lookup longitude")
      sp[driver].add_argument('-mem', type=float, help="Maximum memory allocation when reading from database (GB)")

   return parser, sp


def get_db(args):
   # Set up database
   if args.command == "sim" and args.seed is not None:
      np.random.seed(args.seed)

   model = get_climate_model(args)

   if args.db is None:
      # Don't use args.t as the segment length, because then you never get to join
      # Don't use args.n as the number of segments, because then you never get to join
      if args.dbtype is None or args.dbtype == "random":
         db = wxgen.database.Random(model=model)
      elif args.dbtype == "lorenz63":
         db = wxgen.database.Lorenz63(10, 50, model=wxgen.climate_model.Zero())
      else:
         wxgen.util.error("Cannot understand --dbtype %s" % args.dbtype)
   else:
      db = wxgen.database.Netcdf(args.db, args.vars, model=model, mem=args.mem)

   db.wavelet_levels = args.wavelet_levels

   if args.debug:
      db.info()

   return db


def get_metric(args, db):
   if args.weights is not None:
      weights = np.array(wxgen.util.parse_numbers(args.weights))
      if args.wavelet_levels is not None:
         """
         If we are using wavelets, then the weights in the metric must be repeated such that each
         wavelet component gets the same weight for a given variable.
         """
         NX, NY = db.get_wavelet_size()
         N = int(NX * NY)
         weights = np.repeat(weights, N)

      if args.metric == "rmsd":
         metric = wxgen.metric.Rmsd(weights)
      elif args.metric == "max":
         metric = wxgen.metric.Max(weights)
      elif args.metric == "mad":
         metric = wxgen.metric.Mad(weights)
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


def get_transform(args):
   transform = None
   if args.transform is not None:
      transform = wxgen.transform.get(args.transform)
   return transform


def get_aggregator(args):
   aggregator = None
   if args.aggregator is not None:
      aggregator = wxgen.aggregator.get(args.aggregator)
   return aggregator


def get_module_names(module):
   """
   Returns a list of strings, one for each class in the module
   """
   return [x[0].lower() for x in module.get_all() if "wxgen." + x[0].lower() != module.__name__]


if __name__ == '__main__':
   main()
