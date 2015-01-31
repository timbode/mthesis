# histogram on screen
import os, sys
from numpy import *
from math import *
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

Y=[]; V=[]
counter=0
for s in xrange(SystemData["stats"]):
	if os.path.exists("sysfolder_"+str(s)+"/data/crashed.dat"): continue
	if not os.path.exists("sysfolder_"+str(s)): continue
	liste=[]
	counter=counter+1
	with open("sysfolder_"+str(s)+"/"+SystemData["_DATA_"]+'/init/'+SystemData["system_type"]+'_'+str(0)+'_init_chunk_'+str(1)+'.dat') as f: # p=0, rep=1
		for line in f:
			strang=line.strip().split()
			liste.append(strang)

	Y.append(float(liste[1][0])) # type must be changed to float - keep forgetting that...

	v=[]
	for i in xrange(2):
		v.append(float(liste[i][1]))
	V.append(v)

#print Y

#bins=40
#bin_width=SystemData["a_Y"]/bins
#binning=[q for q in arange(0, SystemData["a_Y"] + bin_width, bin_width)]

# histogram
fig=plt.figure()
plt.hist(Y, bins=100, normed=False, orientation="horizontal", label=str(counter)+" out of "+str(SystemData["stats"]))
plt.title('Location probability on screen')
plt.legend()
plt.xlabel('#')
plt.ylabel('y')
plt.grid(True)
plt.ylim(0.0, SystemData["a_Y"])
#ax=plt.gca()
#ax.invert_xaxis()
#ax.yaxis.set_label_position("right")
#ax.yaxis.tick_right()
fig.savefig('plots/screen_hist.png')
