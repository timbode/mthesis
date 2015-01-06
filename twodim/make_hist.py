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
#---------------------------------------------------------------------------------

system_data.read()
SystemData=system_data.Dict

print "\n========================================================"
print SystemData
print "========================================================\n"

bins_x=system_data.bins_x; bins_y=system_data.bins_y

if  SystemData["system_type"]=='particle':
	# binning
	bin_size_x=float(SystemData["N_X"]-1)/bins_x
	bin_size_y=float(SystemData["N_Y"]-1)/bins_y
	xbins=arange(0, (SystemData["N_X"]-1)+bin_size_x, bin_size_x)
	ybins=arange(0, (SystemData["N_Y"]-1)+bin_size_y, bin_size_y)
elif SystemData["system_type"]=='droplet':
	# binning
	xbins=linspace(0, 1, bins_x+1)
	ybins=linspace(0, 1, bins_y+1)

# 2D histogram
all_counts=zeros((bins_x, bins_y))
for root, _, files in os.walk(SystemData["_DATA_"]+'/hist'):
	for file in files:
		X=[]; Y=[];
		string=file.strip().split('_')
		p=string[2]; rep=string[4][:-4]
		with open(SystemData["_DATA_"]+'/hist/hist_counts_'+p+'_chunk_'+rep+'.dat') as g:
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
if SystemData["N_X"] > SystemData["N_Y"]:
	size_x=size
	size_y=size*(SystemData["N_Y"]-1)/(SystemData["N_X"]-1)
else:
	size_x=size*(SystemData["N_X"]-1)/(SystemData["N_Y"]-1)
	size_y=size
ax=fig.add_axes((0.15, 0.2, size_x, size_y))
plt.title('Location probability')
mesh=ax.pcolormesh(xbins, ybins, all_counts)
fig.colorbar(mesh)
font_size=8
x_text=0.1
fig.text(x_text, 0.04, 'Verlet: '+'dt='+str(SystemData["dt"])+', '+'T='+str(SystemData["T"]), fontsize=font_size)
if  SystemData["system_type"]=='particle':
	fig.text(x_text, 0.07, 'Grid: '+'m='+str(SystemData["m"])+', '+'d='+str(SystemData["d"])+', '+'k='+str(SystemData["k"]), fontsize=font_size)
	fig.text(x_text, 0.1, str.capitalize(SystemData["system_type"])+': '+'M='+str(SystemData["M"])+', '+'D='+str(SystemData["D"]), fontsize=font_size)
	plt.xlim(0, SystemData["N_X"]-1)
	plt.ylim(0, SystemData["N_Y"]-1)
elif SystemData["system_type"]=='droplet':
	fig.text(x_text, 0.07, 'Grid: '+'m='+str(SystemData["m"])+', '+'k='+str(SystemData["k"]), fontsize=font_size)
	fig.text(x_text, 0.1, str.capitalize(SystemData["system_type"])+': '+'M='+str(SystemData["M"])+', '+'f='+str(SystemData["f"]), fontsize=font_size)
	plt.xlim(0, 1)
	plt.ylim(0, 1)
plt.xlabel('x')
plt.ylabel('y')
#fig.tight_layout()
fig.savefig(SystemData["_DATA_"]+'/plots/histogram2D.png')
