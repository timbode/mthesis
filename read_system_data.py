from numpy import *
from math import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

constants=[]
particle_data=[]
lattice_positions=[]
lattice_velocities=[]

with open('particle_data.txt') as f:
	for line in f:
		if line.startswith("#"):
			strang=line.strip().split()
			constants.append(strang)
		else:	
			string=line.strip().split()
			Dict={}
			Dict['b']=int(string[0])
			Dict['pos']=float(string[1])
			Dict['index']=int(string[2])
			Dict['v']=float(string[3])
			Dict['xdot['+string[2]+']']=float(string[4])
			particle_data.append(Dict)
	f.close()
	
constants=constants[:-2]
print constants

'''
#---------------------------------------------------------------------------------

with open('lattice_positions.txt') as f:
	for line in f:
		string=line.strip().split()
		liste=[]
		for s in string:
			liste.append(float(s))
		lattice_positions.append(liste)
	f.close()

#---------------------------------------------------------------------------------

with open('lattice_velocities.txt') as f:
	for line in f:
		string=line.strip().split()
		liste=[]
		for s in string:
			liste.append(float(s))
		lattice_velocities.append(liste)
	f.close()
	
#---------------------------------------------------------------------------------	
'''

def AddInterval(bin_width,c,d,liste,weight,weights): # Interval [c,d]
	for k in xrange(int(floor(c/bin_width)),int(ceil(d/bin_width))):
		liste.append(k*bin_width + bin_width/2)
		weights.append(weight)
	return 0

#---------------------------------------------------------------------------------

L=float(constants[4][2])
N=int(constants[0][2])
steps=int(constants[5][2])
resol=int(constants[6][2])
bins=N/2
bin_width=(N+1)*L/bins
burn_in=steps/5

mean_velo=[]
for i in xrange(0,burn_in):
	#print particle_data[len(particle_data)-i-1]['v']
	mean_velo.append(particle_data[len(particle_data)-i-1]['v'])

print '-------------------------'
print max(mean_velo)
print min(mean_velo)
print sum([abs(q) for q in mean_velo])/burn_in

loc_prob=[]
weights=[]
for j in xrange(0,resol):
	for i in xrange(0,steps - burn_in):
		b=particle_data[j*steps + burn_in + i]['b']
		initial=particle_data[j*steps + burn_in + i-1]['pos']
		final=particle_data[j*steps + burn_in + i]['pos']
		weight=1/abs(particle_data[j*steps + burn_in + i]['v'])
		if b < 0:
			if (b % 2) == 0: # gerade
				AddInterval(bin_width,0,initial,loc_prob,weight,weights)
				AddInterval(bin_width,final,(N+1)*L,loc_prob,weight,weights)
				for k in xrange(0,abs(b)-1): # ganze Strecken
					AddInterval(bin_width,0,(N+1)*L,loc_prob,weight,weights)
			else: # ungerade
				AddInterval(bin_width,0,initial,loc_prob,weight,weights)
				AddInterval(bin_width,0,final,loc_prob,weight,weights)
				for k in xrange(0,abs(b)-1): # ganze Strecken
					AddInterval(bin_width,0,(N+1)*L,loc_prob,weight,weights)
		elif b > 0:
			if (b % 2) == 0: # gerade
				AddInterval(bin_width,initial,(N+1)*L,loc_prob,weight,weights)
				AddInterval(bin_width,0,final,loc_prob,weight,weights)
				for k in xrange(0,abs(b)-1): # ganze Strecken
					AddInterval(bin_width,0,(N+1)*L,loc_prob,weight,weights)
			else: # ungerade
				AddInterval(bin_width,initial,(N+1)*L,loc_prob,weight,weights)
				AddInterval(bin_width,final,(N+1)*L,loc_prob,weight,weights)
				for k in xrange(0,abs(b)-1): # ganze Strecken
					AddInterval(bin_width,0,(N+1)*L,loc_prob,weight,weights)
		elif b == 0:
			if initial < final:
				AddInterval(bin_width,initial,final,loc_prob,weight,weights)
			else:
				AddInterval(bin_width,final,initial,loc_prob,weight,weights)

pp=PdfPages('output.pdf')		
plt.hist(loc_prob, bins=[q for q in arange(0,(N+1)*L + bin_width,bin_width)], normed=0, weights=weights, facecolor='green')
plt.savefig(pp,format='pdf')
plt.clf()
pp.close()
