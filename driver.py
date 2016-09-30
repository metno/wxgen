import sys
import matplotlib.pylab as mpl
from database import *
from trajectory import *
from generator import *

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

database = Database("database.nc")

generator = Generator(database)
#T = 50 # How long is the trajectory?
#N = 10 # How many trajectories?
initial_state = {"T": 0+np.zeros(1, float)}
trajectory = [generator.get(T, initial_state)["T"] for i in range(0, N)]

for tr in trajectory:
   mpl.plot(tr, 'k.-', lw=0.5)
   #mpl.plot(np.diff(tr))
mpl.show()
