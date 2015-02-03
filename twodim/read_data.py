# data analysis for 2D
import os, sys
from numpy import *
from math import *
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams["agg.path.chunksize"]=20000

import system_data

#---------------------------------------------------------------------------------

def Read(file):
	X=[]; Y=[];
	string=file.strip().split('_')
	p=string[1]; rep=string[3][:-4]
	with open(SystemData["_DATA_"]+'/chunks/'+SystemData["system_type"]+'_'+p+'_chunk_'+rep+'.dat') as f:
		for k, line in enumerate(f):
			strang=line.strip().split()
			x=float(strang[0])
			y=float(strang[1])
			X.append(x)
			Y.append(y)
			if float(p) == 0:
				if k % 1000 == 0:
					E.append(float(strang[2]))
					E_grid.append(float(strang[3]))
					E_tot.append(float(strang[2]) + float(strang[3]))

		f.close()

	counts=plt.hist2d(X, Y, bins=[xbins,ybins])[0]
	with open(SystemData["_DATA_"]+'/hist/hist_counts_'+p+'_chunk_'+rep+'.dat', 'w') as g:
		for bracket in counts:
			for element in bracket:
				g.write(str(element)+'\t')
			g.write('\n')

#---------------------------------------------------------------------------------

system_data.read()
SystemData=system_data.Dict

bins_x=system_data.bins_x; bins_y=system_data.bins_y

if  SystemData["system_type"]=='particle':
	# binning
	bin_size_x=float(SystemData["N_X"]-1)/bins_x
	bin_size_y=float(SystemData["N_Y"]-1)/bins_y
	xbins=arange(0, (SystemData["N_X"]-1)+bin_size_x, bin_size_x)
	ybins=arange(0, (SystemData["N_Y"]-1)+bin_size_y, bin_size_y)
elif SystemData["system_type"]=='droplet':
	# binning
	xbins=linspace(0, SystemData["a_X"], bins_x+1) # use number of grid points to derive box size
	ybins=linspace(0, SystemData["a_Y"], bins_y+1)

E=[]; E_grid=[]; E_tot=[];
for root, _, files in os.walk(SystemData["_DATA_"]+'/chunks'):
	for file in files:
		X=[]; Y=[];
		string=file.strip().split('_')
		p=string[1] # p is an artifact
		rep=string[3][:-4]

		with open(SystemData["_DATA_"]+'/chunks/'+SystemData["system_type"]+'_'+p+'_chunk_'+rep+'.dat') as f:
			for k, line in enumerate(f):
				strang=line.strip().split()
				x=float(strang[0])
				y=float(strang[1])
				if (x==0.0 and y==0.0): continue
				if int(p) == 0:
					if k % 1000 == 0:
						E.append(float(strang[2]))
						E_grid.append(float(strang[3]))
						E_tot.append(float(strang[2]) + float(strang[3]))

				# burn-in
				#if int(rep) == 0: continue # watch the types: rep and p have type str

				X.append(x)
				Y.append(y)

			f.close()

		counts=plt.hist2d(X, Y, bins=[xbins,ybins])[0]
		with open(SystemData["_DATA_"]+'/hist/hist_counts_'+p+'_chunk_'+rep+'.dat', 'w') as g:
			for bracket in counts:
				for element in bracket:
					g.write(str(element)+'\t')
				g.write('\n')

#energies
# grid energy
fig=plt.figure(figsize=(20,10))
t_axis=arange(0, len(E))
ax1=plt.subplot(211)
ax1.plot(t_axis, E_grid)
ax1.plot(t_axis, E_tot)
plt.title('Grid energy')
plt.xlabel('t')
plt.ylabel('E')

# droplet energy
ax2=plt.subplot(212)
ax2.plot(t_axis, E)
plt.title('Droplet energy')
plt.xlabel('t')
plt.ylabel('E')
fig.tight_layout()
fig.savefig('data/plots/energies.png')
