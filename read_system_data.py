from numpy import *
from math import *
import matplotlib
import random

#matplotlib.use('Agg')

import matplotlib.pyplot as plt

matplotlib.rcParams["agg.path.chunksize"]=20000
from mpl_toolkits.mplot3d import Axes3D

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

def AddInterval(bin_width,c,d,liste,weight,weights): # Interval [c,d]
	for k in xrange(int(floor(c/bin_width)),int(ceil(d/bin_width))):
		liste.append(k*bin_width + bin_width/2)
		weights.append(weight)
	return 0
	
	
def TimePlots(steps, axis_size, From, To, what='positions'):
	strangs={'positions': particle_positions, 'velocities': particle_velocities, 'lattice_velocities': lattice_velocities}
	for strang in strangs:
			if strang == what:
				liste=strangs[strang]
				
	fig=plt.figure(figsize=(20,20))
	for k in xrange(0,int(steps/axis_size)):
		if k < From/axis_size: continue
		if k >= To/axis_size: break
		time_start=int(k*axis_size)
		time_end=int((k+1)*axis_size)
		
		y_axis=liste[time_start:time_end-1]	
		t_axis=xrange(0, len(y_axis))
		
		ax=plt.subplot(int(To/axis_size) - int(From/axis_size), 1, k+1 - From/axis_size)
		if what == 'positions':
			ax.set_ylim(0,1)
		ax.grid(True)
		ax.plot(t_axis,y_axis)
		xticklabs=ax.get_xticks().tolist()
		xticklabs=[int(time_start + q) for q in xticklabs]
		ax.set_xticklabels(xticklabs)
	fig.savefig('data/'+what+'.png')
	#fig.savefig('data/output_'+str(N)+'_'+str(steps)+'_'+str(delta_t)+'.png')

def ScatterPlot(where, points, what='velocities'):
	strangs={'positions': particle_positions, 'velocities': particle_velocities, 'lattice_velocities': lattice_velocities}
	for strang in strangs:
			if strang == what:
				liste=strangs[strang]
	liste1=liste[:-1]
	liste2=liste[1:]
	liste1=liste1[where:where+points]
	liste2=liste2[where:where+points]
	
	fig=plt.figure()
	plt.scatter(liste1, liste2, s=0.1)
	plt.title('Recursion plot - ' + what)
	if what == 'positions':
		plt.xlim(0,1); plt.ylim(0,1)
	plt.grid(True)
	fig.savefig('data/recur_'+what+'.png')
	
def SmartPlots(steps, axis_size, From, To, M, m):
	liste=[(M*particle_velocities[k] + m*lattice_velocities[k])/(M+m) for k in xrange(0,len(particle_velocities))]
	fig=plt.figure(figsize=(20,200))
	for k in xrange(0,int(steps/axis_size)):
		if k < From/axis_size: continue
		if k >= To/axis_size: break
		time_start=int(k*axis_size)
		time_end=int((k+1)*axis_size)
		
		y_axis=liste[time_start:time_end-1]	
		t_axis=xrange(0, len(y_axis))
		
		ax=plt.subplot(int(To/axis_size) - int(From/axis_size), 1, k+1 - From/axis_size)
		ticker=axis_size/(len(ax.get_xticks()) - 1)
		ax.set_xticklabels(arange(time_start,time_end + ticker, ticker))
		ax.grid(True)
		ax.plot(t_axis,y_axis)
	fig.savefig('data/'+'c-o-m'+'.png')
	
def ThreeD_Plot(From, To):
	fig=plt.figure(figsize=(20,20))
	ax=fig.gca(projection='3d')
	x_axis=particle_velocities[From:To]
	y_axis=lattice_velocities[From:To]
	z_axis=xrange(0,len(x_axis))
	ax.plot(x_axis,y_axis,z_axis,'b--')	
	plt.show()#fig.savefig('data/'+'3D_Plot'+'.png')
		
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

constants=[]
with open('data/system_info.txt') as f:
	for line in f:
		if line.startswith('#'):
			strang=line.strip().split()
			constants.append(strang)
	f.close()
		
N=int(constants[0][2])
M=float(constants[1][2])
m=float(constants[2][2])
L=float(constants[4][2])
steps=int(constants[5][2])
chunk_size=int(constants[6][2])
resol=int(constants[7][2])
delta_t=float(constants[8][2])
print constants

#---------------------------------------------------------------------------------

bins=100
bin_width=(N+1)*L/bins
binning=[q for q in arange(0,(N+1)*L + bin_width,bin_width)]

particle_positions=[]
particle_velocities=[]
lattice_velocities=[]

n1=zeros(bins); n2=zeros(bins); n3=zeros(bins)
chunks_omitted=6
for i in xrange(0,resol):
	initial=0
	for j in xrange(0,steps/chunk_size):
		loc_prob1=[0]; weights1=[0]
		loc_prob2=[0]; weights2=[0]
		loc_prob3=[0]; weights3=[0]
		with open('data/particle_'+str(i)+'_data_'+str(j)+'.txt') as f:
			for k, line in enumerate(f):
				string=line.strip().split()
				particle_positions.append(float(string[1])) # 1 for real positions, 2 for indices
				particle_velocities.append(float(string[3]))
				lattice_velocities.append(float(string[5]))
				
				if j < chunks_omitted:
					continue
			
				b=int(string[0])
				final=float(string[1])
				v=float(string[3])
			
				weight=1/abs(v) # for the moving system: 1/abs(v) --- for the bouncing system: L/abs(v)
			
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
	
		if j >= chunks_omitted:
			n1+=plt.hist(loc_prob1, bins=binning, normed=False, weights=weights1)[0]
			n2+=plt.hist(loc_prob2, bins=binning, normed=False, weights=weights2)[0]
			n3+=plt.hist(loc_prob3, bins=binning, normed=False, weights=weights3)[0]


# open and show 3D plot first
#ThreeD_Plot(300000,310000)
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------	
# histogram
binning=binning[:-1]

fig1=plt.figure()
plt.bar(binning,n2,bin_width)
plt.bar(binning,n3,bin_width, bottom=n2)
plt.bar(binning,n1,bin_width, bottom=n2+n3)
plt.title('Location probability')
plt.xlabel('pos')
plt.ylabel('#')
plt.grid(True)
plt.xlim(0.0, 1.0)
#plt.ylim(0.0, 0.5e8)
fig1.savefig('data/histogram.png')
#fig1.savefig('data/histogram_'+str(N)+'_'+str(steps)+'_'+str(delta_t)+'.png')
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
From=0
To=steps
TimePlots(steps, steps/10, From, To)
TimePlots(steps, steps/10, From, To, what='velocities')
TimePlots(steps, steps/10, From, To, what='lattice_velocities')
#where=random.randint(steps/2,steps)
where=380000#418000# 425000
points=10000
print "where: ", where
ScatterPlot(where, points)
ScatterPlot(where, points, what='positions')
ScatterPlot(where, points, what='lattice_velocities')

#SmartPlots(steps, steps/500, 300000, 400000, M, m)
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
