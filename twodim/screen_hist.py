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
for s in xrange(SystemData["stats"]):
	if not os.path.exists("sysfolder_"+str(s)): continue
	liste=[]
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

bins=20
bin_width=1.0/bins
binning=[q for q in arange(0, 1.0 + bin_width, bin_width)]

# histogram
fig=plt.figure()
plt.hist(Y, bins=binning, normed=False)
plt.title('Location probability on screen')
plt.xlabel('y')
plt.ylabel('#')
plt.grid(True)
plt.xlim(0.0, 1.0)
fig.savefig('plots/screen_hist.png')
