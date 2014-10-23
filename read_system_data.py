from numpy import *
from math import *
import matplotlib
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

#---------------------------------------------------------------------------------

constants=[]
'''
L=float(constants[4][2])
N=int(constants[0][2])
steps=int(constants[5][2])
resol=int(constants[6][2])
'''
N=0
L=0
bin_width=0

bins=50

burn_in=100000

loc_prob1=[]
weights1=[]
loc_prob2=[]
weights2=[]
loc_prob3=[]
weights3=[]

particle_positions=[]
particle_velocities=[]
lattice_velocities=[]
with open('particle_data.txt') as f:
	counter=0
	initial=0
	for k, line in enumerate(f):
		if line.startswith("#"):
			strang=line.strip().split()
			constants.append(strang)
		else:	
			string=line.strip().split()
			particle_positions.append(float(string[1]))
			particle_velocities.append(float(string[3]))
			lattice_velocities.append(float(string[4]))
		
			if k < (counter*int(constants[5][2]) + burn_in + 9):
				continue
				
			if ((k-9) % int(constants[5][2])) == 0:
				counter+=1
				continue
			#print k
				
			b=int(string[0])
			final=float(string[1])
			v=float(string[3])
			
			weight=1/abs(v)
			
			if initial == 0:
				initial=final
				continue
				
			#---------------------------------------------------------------------------------
			N=int(constants[0][2])
			L=float(constants[4][2])
			bin_width=(N+1)*L/bins
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
	
constants=constants[:-2]
print constants

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

pp=PdfPages('output.pdf')		

#---------------------------------------------------------------------------------
# histogram
fig1=plt.figure()
plt.hist([loc_prob1,loc_prob2,loc_prob3], bins=[q for q in arange(0,(N+1)*L + bin_width,bin_width)], normed=False, weights=[weights1,weights2,weights3], stacked=True)
plt.title('Location probability')
plt.xlabel('pos')
plt.ylabel('#')
plt.grid(True)
#plt.ylim(0.0, 0.7e+5)
pp.savefig(fig1)
#plt.clf()

time_plots-0
if time_plots==1:
	#---------------------------------------------------------------------------------
	# positions
	y_axis=particle_positions		
	t_axis=xrange(0, len(y_axis))

	fig2=plt.figure(figsize=(20,5))
	plt.title('Trajectory')
	plt.xlabel('t')
	plt.ylabel('pos')
	plt.plot(t_axis,y_axis)
	plt.grid(True)
	#plt.axis([0,500000,0,1])
	pp.savefig(fig2)

	#---------------------------------------------------------------------------------

	# velocities
	y_axis=particle_velocities		
	t_axis=xrange(0, len(y_axis))

	fig3=plt.figure(figsize=(20,5))
	plt.title('Velocity (post)')
	plt.xlabel('t')
	plt.ylabel('v')
	plt.plot(t_axis,y_axis)
	plt.grid(True)
	#plt.axis([0,500000,0,1])
	pp.savefig(fig3)

	#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
recursion_plots=1
start=1100000
end=1150000

if recursion_plots==1:
	#---------------------------------------------------------------------------------
	# recursion plot positions
	liste1=particle_positions[:-1]
	liste2=particle_positions[1:]
	liste1=liste1[start:end]
	liste2=liste2[start:end]

	fig4=plt.figure()
	plt.scatter(liste1,liste2,s=0.1)
	plt.title('Recursion plot - positions of collisions')
	plt.xlabel('pos_i')
	plt.ylabel('pos_i+1')
	plt.axis([0,1,0,1])
	plt.grid(True)
	pp.savefig(fig4)
	#---------------------------------------------------------------------------------
	# recursion plot velocities
	liste1=particle_velocities[:-1]
	liste2=particle_velocities[1:]
	liste1=liste1[start:end]
	liste2=liste2[start:end]

	fig5=plt.figure()
	plt.title('Recursion plot - velocities (post)')
	plt.xlabel('v_i')
	plt.ylabel('v_i+1')
	plt.scatter(liste1,liste2,s=0.1)
	plt.grid(True)
	pp.savefig(fig5)
	#---------------------------------------------------------------------------------
	# recursion plot lattice velocities
	liste1=lattice_velocities[:-1]
	liste2=lattice_velocities[1:]
	liste1=liste1[start:end]
	liste2=liste2[start:end]

	fig6=plt.figure()
	plt.scatter(liste1,liste2,s=0.1)
	plt.title('Recursion plot - velocity of current lattice mass (post)')
	plt.xlabel('xdot_i')
	plt.ylabel('xdot_i+1')
	plt.grid(True)
	pp.savefig(fig6)
#---------------------------------------------------------------------------------

pp.close()

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
