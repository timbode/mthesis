from numpy import *
from math import *
import matplotlib
import random
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams["agg.path.chunksize"]=20000

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

def AddInterval(bin_width,c,d,liste,weight,weights): # Interval [c,d]
	for k in xrange(int(floor(c/bin_width)),int(ceil(d/bin_width))):
		liste.append(k*bin_width + bin_width/2)
		weights.append(weight)
	return 0
	
	
def TimePlots(steps, axis_size, what='positions'):
	strangs={'positions': particle_positions, 'velocities': particle_velocities, 'lattice_velocities': lattice_velocities}
	for strang in strangs:
			if strang == what:
				liste=strangs[strang]
				
	fig=plt.figure(figsize=(20,20))
	for k in xrange(0,int(steps/axis_size)):
		time_start=int(k*axis_size)
		time_end=int((k+1)*axis_size)
		
		y_axis=liste[time_start:time_end-1]	
		t_axis=xrange(0, len(y_axis))
		
		ax=plt.subplot(int(steps/axis_size),1,k+1)
		ticker=axis_size/(len(ax.get_xticks()) - 1)
		ax.set_xticklabels(arange(time_start,time_end + ticker, ticker))
		if what == 'positions':
			ax.set_ylim(0,1)
		ax.grid(True)
		ax.plot(t_axis,y_axis)
	fig.savefig('data/'+what+'.png')
	#fig.savefig('data/output_'+str(N)+'_'+str(steps)+'_'+str(delta_t)+'.png')

def ScatterPlots(where, points, what='velocities'):
	strangs={'positions': particle_positions, 'velocities': particle_velocities, 'lattice_velocities': lattice_velocities}
	for strang in strangs:
			if strang == what:
				liste=strangs[strang]
	liste1=liste[:-1]
	liste2=liste[1:]
	print where, points
	liste1=liste1[where:where+points]
	liste2=liste2[where:where+points]
	
	fig=plt.figure()
	plt.scatter(liste1, liste2, s=0.1)
	plt.title('Recursion plot - ' + what)
	#plt.xlim(0,1); plt.ylim(0,1)
	plt.grid(True)
	fig.savefig('data/recur_'+what+'.png')

#---------------------------------------------------------------------------------

constants=[]
with open('data/system_info.txt') as f:
	for line in f:
		if line.startswith('#'):
			strang=line.strip().split()
			constants.append(strang)
	f.close()
		
N=int(constants[0][2])
L=float(constants[4][2])
steps=int(constants[5][2])
resol=int(constants[6][2])
delta_t=float(constants[7][2])
print constants

bins=100
bin_width=(N+1)*L/bins
burn_in=1000

loc_prob1=[]; weights1=[]
loc_prob2=[]; weights2=[]
loc_prob3=[]; weights3=[]

particle_positions=[]
particle_velocities=[]
lattice_velocities=[]
for i in xrange(0,resol):
	with open('data/particle_data_'+str(i)+'.txt') as f:
		initial=0
		for j, line in enumerate(f):
			string=line.strip().split()
			particle_positions.append(float(string[1])) # 1 for real positions, 2 for indices
			particle_velocities.append(float(string[3]))
			lattice_velocities.append(float(string[5]))
		
			if j < burn_in:
				continue
			
			b=int(string[0])
			final=float(string[1])
			v=float(string[3])
			
			weight=1/abs(v)
			
			if initial == 0:
				initial=final
				continue
				
			#---------------------------------------------------------------------------------
			if b < 0:
				if (b % 2) != 0: # ungerade
					AddInterval(bin_width,0,initial,loc_prob2,weight,weights2)
					AddInterval(bin_width,0,final,loc_prob2,weight,weights2)
					for k in xrange(0,abs(b)-1): # ganze Strecken
						AddInterval(bin_width,0,(N+1)*L,loc_prob2,weight,weights2)
			
				else: # gerade
					AddInterval(bin_width,0,initial,loc_prob2,weight,weights2)
					AddInterval(bin_width,final,(N+1)*L,loc_prob2,weight,weights2)
					for k in xrange(0,abs(b)-1): # ganze Strecken
						AddInterval(bin_width,0,(N+1)*L,loc_prob2,weight,weights2)
			
			if b > 0:
				if (b % 2) != 0: # ungerade
					AddInterval(bin_width,initial,(N+1)*L,loc_prob3,weight,weights3)
					AddInterval(bin_width,final,(N+1)*L,loc_prob3,weight,weights3)
					for k in xrange(0,abs(b)-1): # ganze Strecken
						AddInterval(bin_width,0,(N+1)*L,loc_prob3,weight,weights3)
			
				else: # gerade
					AddInterval(bin_width,initial,(N+1)*L,loc_prob3,weight,weights3)
					AddInterval(bin_width,0,final,loc_prob3,weight,weights3)
					for k in xrange(0,abs(b)-1): # ganze Strecken
						AddInterval(bin_width,0,(N+1)*L,loc_prob3,weight,weights3)
			
			if b == 0:
				if initial < final:
					AddInterval(bin_width,initial,final,loc_prob1,weight,weights1)
				else:
					AddInterval(bin_width,final,initial,loc_prob1,weight,weights1)
			#---------------------------------------------------------------------------------
			
			# Nicht vergessen:
			initial=final
	f.close()
	
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

	

#---------------------------------------------------------------------------------
# histogram
#pp=PdfPages('data/histogram.pdf')	
fig1=plt.figure()
plt.hist([loc_prob1,loc_prob2,loc_prob3], bins=[q for q in arange(0,(N+1)*L + bin_width,bin_width)], normed=False, weights=[weights1,weights2,weights3],stacked=True)
plt.title('Location probability')
plt.xlabel('pos')
plt.ylabel('#')
plt.grid(True)
plt.xlim(0.0, 1.0)
fig1.savefig('data/histogram.png')
#fig1.savefig('data/histogram_'+str(N)+'_'+str(steps)+'_'+str(delta_t)+'.png')
#plt.clf()
#pp.close()
#---------------------------------------------------------------------------------


TimePlots(steps, steps/5, what='velocities')
ScatterPlots(random.randint(0,steps), 10000)
ScatterPlots(0, 10000, what='positions')
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
