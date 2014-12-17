# data analysis for 2D
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

dim=int(constants[0][2])
N_X=int(constants[1][2])
N_Y=int(constants[2][2])
repeat=int(constants[10][2])
stats=int(constants[12][2])

# binning
xbins=arange(0, N_X)
ybins=arange(0, N_Y)

for p in xrange(0, stats):
	if p==0:
		E=[]; E_grid=[]; E_tot=[];
	for rep in xrange(0, repeat):
		X=[]; Y=[];
		with open('data/chunks/particle_'+str(p)+'_chunk_'+str(rep)+'.dat') as f:
			for line in f:
				strang=line.strip().split()
				X.append(float(strang[0]))
				Y.append(float(strang[1]))
				if p==0:
					E.append(float(strang[2]))
					E_grid.append(float(strang[3]))
					E_tot.append(float(strang[2]) + float(strang[3]))
			
			f.close()
			
		counts=plt.hist2d(X, Y, bins=[xbins,ybins])[0]	
		with open('data/hist/hist_counts_'+str(p)+'_chunk_'+str(rep)+'.dat', 'w') as g:
			for bracket in counts:
				for element in bracket:
					g.write(str(element)+'\t')
				g.write('\n')

# energy
fig=plt.figure(figsize=(20,10))
plt.title('Energies')
t_axis=arange(0, len(E))
#plt.plot(t_axis, E)
plt.plot(t_axis, E_grid)
plt.plot(t_axis, E_tot)
plt.xlabel('t')
plt.ylabel('E')
fig.savefig('data/plots/energies.png')
