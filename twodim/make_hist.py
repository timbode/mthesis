# make historgram
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
repeat=int(constants[10][2])
stats=int(constants[12][2])


# 2D histogram
all_counts=zeros((N_Y-1, N_X-1))
for p in xrange(0, stats):
	for rep in xrange(0, repeat):
		with open('data/hist/hist_counts_'+str(p)+'_chunk_'+str(rep)+'.dat') as f:
			counts=[]
			for line in f:
				strang=line.strip().split()
				strang=[float(q) for q in strang]
				counts.append(strang)
			counts=transpose(array(counts))
			all_counts+=counts
			
# binning
xbins=arange(0, N_X)
ybins=arange(0, N_Y)
xbins, ybins=meshgrid(xbins, ybins)

fig=plt.figure()
plt.title('Location probability')
plt.pcolormesh(xbins, ybins, all_counts)
plt.colorbar()
plt.xlim(0,N_X-1)
plt.ylim(0,N_Y-1)
plt.xlabel('x')
plt.ylabel('y')
fig.savefig('data/plots/histogram2D.png')
