# data analysis
from numpy import *
from math import *
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

#matplotlib.rcParams["agg.path.chunksize"]=20000

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

constants=[]
with open('data/system.dat') as f:
	for line in f:
		if line.startswith('#'):
			strang=line.strip().split()
			constants.append(strang)
	f.close()

print "================================================="
print constants
print "================================================="


N_X=int(constants[1][2])
N_Y=int(constants[2][2])
stats=int(constants[11][2])

for p in xrange(0, stats):
	X=[]; Y=[]; Z=[];
	with open('data/particle_'+str(p)+'.dat') as f:
		for line in f:
			strang=line.strip().split()
			strang.append("0.0")
			X.append(float(strang[0]))
			Y.append(float(strang[1]))
			Z.append(float(strang[2]))			
			
		f.close()
		
# 2D histogram
fig=plt.figure()
plt.title('Location probability')
plt.hist2d(X, Y, bins=20)
plt.colorbar()
plt.xlim(0,N_X-1);
plt.ylim(0,N_Y-1)
plt.xlabel('x')
plt.ylabel('y')
fig.savefig('data/histogram2D.png')
