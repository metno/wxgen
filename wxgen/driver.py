import sys
import matplotlib.pylab as mpl
import database
from trajectory import *
from generator import *
import verif

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
      V = 1
      db = database.Random(N, T, V)
   else:
      db_filename = sys.argv[3]
      db = database.Netcdf(db_filename)

   var = db.vars()[0]

   if 0:
      for i in range(0, db.size()):
         mpl.plot(db.get(i)["air_temperature_2m"], 'k.-', lw=0.5)
      mpl.show()
      sys.exit()

   generator = Generator(db)
   initial_state = np.array([275, 0, 0, 0])
   initial_state = np.array([275])
   trajectories = generator.get(N, T, initial_state)

   v = 0
   for tr in trajectories:
      mpl.plot(tr[:,v], 'k.-', lw=0.5)
      #mpl.plot(np.diff(tr))
      #print tr[var]

   v = verif.Change()
   #mpl.plot(np.linspace(0.5, T-0.5, T-1), v.compute(trajectories), 'k.-')
   mpl.show()

if __name__ == '__main__':
       main()
