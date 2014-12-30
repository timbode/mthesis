# data analysis for 2D
import os
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

dim=int(constants[0][2]) # could also iterate a dictionary
N_X=int(constants[1][2])
N_Y=int(constants[2][2])
k=float(constants[4][2])
m=float(constants[5][2])
d=float(constants[6][2])
M=float(constants[7][2])
D=float(constants[8][2])
steps=int(constants[9][2])
repeat=int(constants[10][2])
dt=float(constants[11][2])
stats=int(constants[12][2])
T=steps*repeat*dt

bins_x=50; bins_y=50

# binning
xbins=linspace(0, N_X, bins_x)
ybins=linspace(0, N_Y, bins_y)

E=[]; E_grid=[]; E_tot=[];
for root, _, files in os.walk('data/chunks'):
	for file in files:
		X=[]; Y=[];
		string=file.strip().split('_')
		p=string[1]; rep=string[3][:-4]
		with open('data/chunks/particle_'+p+'_chunk_'+rep+'.dat') as f:
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
'''
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
'''
