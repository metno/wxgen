import sys
import matplotlib.pylab as mpl
import database
from trajectory import *
from generator import *
import verif

if len(sys.argv) < 3:
   print "Weather generator"
   print "usage: wxgen N T"
   print ""
   print "Arguments:"
   print "  N: Number of trajectories"
   print "  T: Number of days in trajectory"

   sys.exit()

N = int(sys.argv[1])
T = int(sys.argv[2])
#np.random.seed(1)

db = database.Random()

generator = Generator(db)
initial_state = {"T": 0+np.zeros(1, float)}
trajectories = [generator.get(T, initial_state) for i in range(0, N)]

#for tr in trajectories:
   #mpl.plot(tr["T"], 'k.-', lw=0.5)
   #mpl.plot(np.diff(tr))

v = verif.Change()
mpl.plot(np.linspace(0.5, T-0.5, T-1), v.compute(trajectories), 'k.-')
mpl.show()
