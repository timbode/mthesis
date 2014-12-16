# data analysis
from numpy import *
from math import *
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

matplotlib.rcParams["agg.path.chunksize"]=20000

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

dim=int(constants[0][2])
N_X=int(constants[1][2])
N_Y=int(constants[2][2])
stats=int(constants[11][2])

for p in xrange(0, stats):
	X=[]; Y=[]; Z=[];
	E=[]; E_grid=[]; E_sum=[];
	with open('data/particle_'+str(p)+'.dat') as f:
		for line in f:
			strang=line.strip().split()
			X.append(float(strang[0]))
			Y.append(float(strang[1]))
			#Z.append(float(strang[2]))			
			E.append(float(strang[2]))
			E_grid.append(float(strang[3]))
			E_sum.append(float(strang[2]) + float(strang[3]))
			
		f.close()
		
# 2D histogram
fig=plt.figure()
plt.title('Location probability')
plt.hist2d(X, Y, bins=100)
plt.colorbar()
plt.xlim(0,N_X-1)
plt.ylim(0,N_Y-1)
plt.xlabel('x')
plt.ylabel('y')
fig.savefig('data/histogram2D.png')

# energy
fig=plt.figure(figsize=(20,10))
plt.title('Energies')
t_axis=arange(0, len(E))
plt.plot(t_axis, E)
plt.plot(t_axis, E_grid)
plt.plot(t_axis, E_sum)
plt.xlabel('t')
plt.ylabel('E')
fig.savefig('data/energies.png')
