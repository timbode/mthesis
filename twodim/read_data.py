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

dim=int(constants[0][2]) # could as well iterate a dictionary
N_X=int(constants[1][2])
N_Y=int(constants[2][2])
k=float(constants[4][2])
m=float(constants[5][2])
d=float(constants[6][2])
L=float(constants[7][2])
M=float(constants[7+1][2])
D=float(constants[8+1][2])
f=float(constants[9+1][2])

steps=int(constants[9+1+1][2])
repeat=int(constants[10+1+1][2])
dt=float(constants[11+1+1][2])
T=steps*repeat*dt

stats=int(constants[12+1+1][2])
system_type=str(constants[13+1+1][2])

bins_x=50; bins_y=50

if  system_type=='particle':
	# binning
	bin_size_x=float(N_X-1)/bins_x
	bin_size_y=float(N_Y-1)/bins_y
	xbins=arange(0, (N_X-1)+bin_size_x, bin_size_x)
	ybins=arange(0, (N_Y-1)+bin_size_y, bin_size_y)
elif system_type=='droplet':
	# binning
	xbins=linspace(0, 1, bins_x+1)
	ybins=linspace(0, 1, bins_y+1)

E=[]; E_grid=[]; E_tot=[];
for root, _, files in os.walk('data/chunks'):
	for file in files:
		X=[]; Y=[];
		string=file.strip().split('_')
		p=string[1]; rep=string[3][:-4]
		with open('data/chunks/'+system_type+'_'+p+'_chunk_'+rep+'.dat') as f:
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
		with open('data/hist/hist_counts_'+p+'_chunk_'+rep+'.dat', 'w') as g:
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
