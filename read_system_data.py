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
			constants.append(line.strip())
		else:	
			string=line.strip().split()
			Dict={}
			Dict['pos']=float(string[0])
			Dict['index']=int(string[1])
			Dict['v']=float(string[2])
			Dict['xdot['+string[1]+']']=float(string[3])
			particle_data.append(Dict)
	f.close()	

with open('lattice_positions.txt') as f:
	for line in f:
		string=line.strip().split()
		liste=[]
		for s in string:
			liste.append(float(s))
		lattice_positions.append(liste)
	f.close()

with open('lattice_velocities.txt') as f:
	for line in f:
		string=line.strip().split()
		liste=[]
		for s in string:
			liste.append(float(s))
		lattice_velocities.append(liste)
	f.close()
	
#------------------------------------------------------------------------------------------------------------------------------------------------------

num=2000
mean_velo=[]
for i in xrange(0,num):
	#print particle_data[len(particle_data)-i-1]['v']
	mean_velo.append(abs(particle_data[len(particle_data)-i-1]['v']))
	
t=arange(0,len(mean_velo))

print constants
print '-------------------------'
#print max(mean_velo)
#print min(mean_velo)
print sum([q for q in mean_velo])/num

#------------------------------------------------------------------------------------------------------------------------------------------------------

time=[]
for i in xrange(1, len(lattice_positions)+1):
	lattice_positions[i-1]=[q+1 + lattice_positions[i-1][q] for q in xrange(0,len(lattice_positions[i-1]))]
	times=[i for q in xrange(0,len(lattice_positions[i-1]))]
	time.append(times)
	
lattice_positions=[q for sublist in lattice_positions for q in sublist]
time=[q for sublist in time for q in sublist]

pp=PdfPages('output.pdf')

plt.plot(t,mean_velo)
plt.savefig(pp,format='pdf')

plt.clf()

plt.scatter(lattice_positions,time,s=2.0)
plt.savefig(pp,format='pdf')

pp.close()
