from numpy import *
from math import *
import matplotlib.pyplot as plt

only_lattice_positions=[]
only_lattice_velocities=[]

with open('only_lattice_positions.txt') as f:
	for line in f:
		string=line.strip().split()
		liste=[]
		for s in string:
			liste.append(float(s))
		only_lattice_positions.append(liste)
	f.close()

#---------------------------------------------------------------------------------

with open('only_lattice_velocities.txt') as f:
	for line in f:
		string=line.strip().split()
		liste=[]
		for s in string:
			liste.append(float(s))
		only_lattice_velocities.append(liste)
	f.close()
	
#---------------------------------------------------------------------------------	

time=[]
for i in xrange(1, len(only_lattice_positions)+1):
	only_lattice_positions[i-1]=[q+1 + only_lattice_positions[i-1][q] for q in xrange(0,len(only_lattice_positions[i-1]))]
	times=[i for q in xrange(0,len(only_lattice_positions[i-1]))]
	time.append(times)
	
only_lattice_positions=[q for sublist in only_lattice_positions for q in sublist]
time=[q for sublist in time for q in sublist]
plt.scatter(only_lattice_positions,time,s=2.0)
plt.show()

