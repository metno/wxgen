import sys
import argparse
import numpy as np
import matplotlib.pylab as mpl
import wxgen.database
import wxgen.output_db
import wxgen.version
try:
   from netCDF4 import Dataset as netcdf
except:
   from scipy.io.netcdf import netcdf_file as netcdf

def run(argv):
   parser = argparse.ArgumentParser(description="Weather generator database")
   parser.add_argument('action', type=str, metavar="action", help="Number of trajectories")
   parser.add_argument('-o', type=str, default=None, help="Output filename", dest="output_filename")
   parser.add_argument('-i', type=str, default=None, help="Input filenames", nargs="+", dest="input_filenames")
   parser.add_argument('--db', type=str, default=None, help="Filename of NetCDF database", required=False)
   parser.add_argument('--type', type=str, default="timeseries", help="Output type (text or plot)")

   args = parser.parse_args()

   if args.action == "plot":
      # Set up database
      db = wxgen.database.Netcdf(args.db)

      output = wxgen.output_db.get(args.type)(db, args.output_filename)
      output.plot()

   elif args.action == "create":
      if args.input_filenames is None:
         wxgen.util.error("Missing input filenames")
      create_database(args.input_filenames, args.output_filename)

def create_database(input_filenames, output_filename, vars=["air_temperature_2m", "x_wind_10m", "y_wind_10m"]):
   """
   Creates a netcdf file
   
   vars     Which variables to output
   """

   # Which gridpoints to averages across
   xrange = range(100, 130)
   yrange = range(82, 107)

   values = dict()
   D = len(input_filenames) # Number of dates
   O = 10 # Number of leadtimes
   E = 51 # Number of ensemble members

   for var in vars:
      values[var] = np.zeros([D, O, E], float)
   dates = np.zeros(D, int)

   # Loop over files
   for d in range(0, D):
      ifile = netcdf(input_filenames[d], 'r')
      dates[d] = ifile.variables["time"][0]
      lats = ifile.variables["latitude"]
      lons = ifile.variables["longitude"]
      for var in vars:
         data = ifile.variables[var][:, 0, :, xrange, yrange]
         temp = np.mean(np.mean(data, axis=3), axis=2)
         for t in range(0, 10):
            I = range(t * 4, (t + 1) * 4)
            values[var][d, t, :] = np.mean(temp[I, :], axis=0)
      ifile.close()

   # Write to output file
   ofile = netcdf(output_filename, 'w')
   ofile.createDimension("date", D)
   ofile.createDimension("leadtime", O)
   ofile.createDimension("member", E)

   vVars = dict()
   vDate = ofile.createVariable("date", "i4", ("date"))
   vDate[:] = dates
   vOffset = ofile.createVariable("leadtime", "i4", ("leadtime"))
   vOffset[:] = range(0, O)
   for var in vars:
      vVars[var] = ofile.createVariable(var, "f4", ("date", "leadtime", "member"))
      vVars[var][:] = values[var]
   ofile.close()

if __name__ == '__main__':
       main()
