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
	#if os.path.exists("sysfolder_"+str(s)+"/data/touched.dat"): continue
	if not os.path.exists("sysfolder_"+str(s)): continue
	liste=[]
	counter=counter+1
	with open("sysfolder_"+str(s)+"/"+SystemData["_DATA_"]+'/init/'+SystemData["system_type"]+'_'+str(0)+'_init_chunk_'+str(1)+'.dat') as f: # p=0, rep=1
		for line in f:
			strang=line.strip().split()
			liste.append(strang)

	Y.append(float(liste[1][0])) # type must be changed to float - keep forgetting that...

	v=[]
	v.append(float(liste[0][1]))
	v.append(float(liste[1][1]))
	norm=v[0]*v[0] + v[1]*v[1]
	norm=sqrt(norm)
	v[0]=v[0]/norm; v[1]=v[1]/norm;

	V.append(v)

#alpha=[(q[1]/abs(q[1]))*acos(q[0]) for q in V]
alpha=[q[1] for q in V]

# histogram
fig=plt.figure()
plt.hist(Y, bins=40, normed=False, orientation="horizontal", label=str(counter)+" out of "+str(SystemData["stats"]))
plt.title('Location probability on screen')
plt.legend()
plt.xlabel('#')
plt.ylabel('y')
#plt.yticks(arange(0, SystemData["a_Y"]+0.1, 0.1))
plt.grid(True)
plt.ylim(0.0, SystemData["a_Y"])
#ax=plt.gca()
#ax.invert_xaxis()
#ax.yaxis.set_label_position("right")
#ax.yaxis.tick_right()
fig.tight_layout()
fig.savefig('plots/screen_hist.png')

'''
# position against velocity angle
fig=plt.figure()
plt.scatter(Y, alpha)
plt.title('Position against velocity angle')
#plt.legend()
plt.grid(True)
plt.xlabel('y')
plt.ylabel('alpha')
plt.xticks(arange(0, SystemData["a_Y"]+0.1, 0.1))
plt.xlim(0.0, SystemData["a_Y"])
#ax=plt.gca()
#ax.invert_xaxis()
#ax.yaxis.set_label_position("right")
#ax.yaxis.tick_right()
fig.tight_layout()
fig.savefig('plots/pos_angle.png')
'''
