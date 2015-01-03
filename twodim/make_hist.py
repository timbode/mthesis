# make historgram
import os
import string
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

# 2D histogram
all_counts=zeros((bins_x, bins_y))
for root, _, files in os.walk('data/hist'):
	for file in files:
		X=[]; Y=[];
		string=file.strip().split('_')
		p=string[2]; rep=string[4][:-4]
		with open('data/hist/hist_counts_'+p+'_chunk_'+rep+'.dat') as g:
			counts=[]
			for line in g:
				strang=line.strip().split()
				strang=[float(q) for q in strang]
				counts.append(strang)
			counts=transpose(array(counts))
			all_counts+=counts


# meshgrid
xbins, ybins=meshgrid(xbins, ybins)
fig=plt.figure()
size=0.7
if N_X > N_Y:
	size_x=size
	size_y=size*(N_Y-1)/(N_X-1)
else:
	size_x=size*(N_X-1)/(N_Y-1)
	size_y=size
ax=fig.add_axes((0.15, 0.2, size_x, size_y))
plt.title('Location probability')
mesh=ax.pcolormesh(xbins, ybins, all_counts)
fig.colorbar(mesh)
font_size=8
x_text=0.1
fig.text(x_text, 0.04, 'Verlet: '+'dt='+str(dt)+', '+'T='+str(T), fontsize=font_size)
if  system_type=='particle':
	fig.text(x_text, 0.07, 'Grid: '+'m='+str(m)+', '+'d='+str(d)+', '+'k='+str(k), fontsize=font_size)
	fig.text(x_text, 0.1, str.capitalize(system_type)+': '+'M='+str(M)+', '+'D='+str(D), fontsize=font_size)
	plt.xlim(0, N_X-1)
	plt.ylim(0, N_Y-1)
elif system_type=='droplet':
	fig.text(x_text, 0.07, 'Grid: '+'m='+str(m)+', '+'k='+str(k), fontsize=font_size)
	fig.text(x_text, 0.1, str.capitalize(system_type)+': '+'M='+str(M)+', '+'f='+str(f), fontsize=font_size)
	plt.xlim(0, 1)
	plt.ylim(0, 1)
plt.xlabel('x')
plt.ylabel('y')
#fig.tight_layout()
fig.savefig('data/plots/histogram2D.png')
