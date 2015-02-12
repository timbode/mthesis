# grid energy

import os, sys
from numpy import *
from math import *
from mpl_toolkits.mplot3d import axes3d
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams["agg.path.chunksize"]=20000

import system_data

#---------------------------------------------------------------------------------
# needs 'data/system.dat' and folder 'plots'
#---------------------------------------------------------------------------------

system_data.read()
SystemData=system_data.Dict

bins_x=system_data.bins_x; bins_y=system_data.bins_y

#---------------------------------------------------------------------------------
def Index(X, Y, Z, Alpha):
	# Watch out: SystemData["N_Z"]=N_[2] omitted
	return X*SystemData["N_Y"] + Y + Z + Alpha*SystemData["N_X"]*SystemData["N_Y"]
#---------------------------------------------------------------------------------

# 2D histogram
for s in xrange(SystemData["stats"]):
	if not os.path.exists("sysfolder_"+str(s)): continue
	rdot=[]
	with open("sysfolder_"+str(s)+"/"+"burn_off"+".dat") as g:
		for line in g:
			strang=line.strip().split()
			v=float(strang[0])
			rdot.append(v)


	# looping...
	E=zeros((SystemData["N_Y"], SystemData["N_X"]))
	for x in xrange(SystemData["N_X"]):
		for y in xrange(SystemData["N_Y"]):
			for z in xrange(1):
				e=0;
				for alpha in xrange(2):
					#e+=0.5*rdot[Index(x, y, z, alpha)]*rdot[Index(x, y, z, alpha)]
					if alpha==0:
						temp=rdot[Index(x, y, z, alpha)]#*rdot[Index(x, y, z, alpha)]
						if temp==0: continue
						if abs(temp) > 0.01: continue
						#if temp > 0: temp=temp+exp(abs(temp))
						#else: temp=temp-exp(abs(temp))
					
						e+=temp
				E[y][x]+=e

	xbins=arange(0, SystemData["N_X"]) # use number of grid points to derive box size
	ybins=arange(0, SystemData["N_Y"])
	X, Y=meshgrid(xbins, ybins)
	E=array(E)

	fig=plt.figure()
	plt.title('Grid kinetic energy')
	plt.pcolormesh(X, Y, E, cmap=plt.cm.Blues)# Blues, jet
	plt.grid(True, color='white')
	plt.xlabel('x')
	plt.ylabel('y')
	fig.tight_layout()
	fig.savefig("sysfolder_"+str(s)+"/"+'data/plots/grid_energy.png')
