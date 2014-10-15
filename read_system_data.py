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
bins=100
bin_width=(N+1)*L/bins
burn_in=100000

loc_prob1=[]
weights1=[]
loc_prob2=[]
weights2=[]
loc_prob3=[]
weights3=[]
for j in xrange(0,resol):
	for i in xrange(0,steps - burn_in):
		b=particle_data[j*steps + burn_in + i]['b']
		initial=particle_data[j*steps + burn_in + i-1]['pos']
		final=particle_data[j*steps + burn_in + i]['pos']
		weight=1/abs(particle_data[j*steps + burn_in + i]['v'])
		'''		
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
			
		'''
		if b > 0:
			if (b % 2) != 0: # ungerade
				AddInterval(bin_width,initial,(N+1)*L,loc_prob3,weight,weights3)
				AddInterval(bin_width,final,(N+1)*L,loc_prob3,weight,weights3)
				print weight, particle_data[j*steps + burn_in + i]['v']
				for k in xrange(0,abs(b)-1): # ganze Strecken
					AddInterval(bin_width,0,(N+1)*L,loc_prob3,weight,weights3)
			
			else: # gerade
				AddInterval(bin_width,initial,(N+1)*L,loc_prob3,weight,weights3)
				AddInterval(bin_width,0,final,loc_prob3,weight,weights3)
				for k in xrange(0,abs(b)-1): # ganze Strecken
					AddInterval(bin_width,0,(N+1)*L,loc_prob3,weight,weights3)
			
		if b == 0:
			print "__________________________", weight, particle_data[j*steps + burn_in + i]['v']
			if initial < final:
				AddInterval(bin_width,initial,final,loc_prob1,weight,weights1)
			else:
				AddInterval(bin_width,final,initial,loc_prob1,weight,weights1)



pp=PdfPages('output.pdf')		
plt.figure()
plt.hist([loc_prob1,loc_prob2,loc_prob3], bins=[q for q in arange(0,(N+1)*L + bin_width,bin_width)], normed=False, weights=[weights1,weights2,weights3])
#plt.ylim(0.0, 0.7e+5)
plt.savefig(pp,format='pdf')

plt.clf()

liste=particle_data[15550:15570]		
t_axis=xrange(0, len(liste))
x_axis=[q['pos'] for q in liste] 

plt.plot(t_axis,x_axis)
plt.savefig(pp, format='pdf')

plt.clf()

liste=particle_data[15570:15670]		
t_axis=xrange(0, len(liste))
x_axis=[q['pos'] for q in liste] 

plt.plot(t_axis,x_axis)
plt.savefig(pp, format='pdf')

plt.clf()

liste=particle_data[15670:19900]		
t_axis=xrange(0, len(liste))
x_axis=[q['pos'] for q in liste] 

plt.plot(t_axis,x_axis)
plt.savefig(pp, format='pdf')

pp.close()
