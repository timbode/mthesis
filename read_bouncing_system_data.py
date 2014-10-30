from numpy import *
from math import *
import matplotlib
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

#---------------------------------------------------------------------------------

constants=[]
with open('data/bouncing_system_info.txt') as f:
	for line in f:
		if line.startswith('#'):
			strang=line.strip().split()
			constants.append(strang)
	f.close()
		
N=int(constants[0][2])
L=float(constants[4][2])
steps=int(constants[5][2])
resol=int(constants[6][2])
print constants

bins=50
bin_width=(N+1)*L/bins
burn_in=10000

loc_prob=[]
weights=[]

particle_positions=[]
particle_velocities=[]
lattice_velocities_pre=[]
for i in xrange(0,resol):
	with open('data/bouncing_particle_data_'+str(i)+'.txt') as f:
		initial=0
		for j, line in enumerate(f):
			string=line.strip().split()
			particle_positions.append(L*(float(string[0])+1))
			particle_velocities.append(float(string[1]))
			lattice_velocities_pre.append(float(string[2]))
			
			if j < burn_in:
				continue		
				
			final=L*(int(string[0])+1)
			v=float(string[1])
			
			weight=L/abs(v)
			
			if initial == 0:
				initial=final
				continue

			if initial < final:
				AddInterval(bin_width,initial,final,loc_prob,weight,weights)
			else:
				AddInterval(bin_width,final,initial,loc_prob,weight,weights)				
			
			# Nicht vergessen:
			initial=final
	f.close()
	
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

pp=PdfPages('data/bouncing_output.pdf')		

#---------------------------------------------------------------------------------
# histogram
fig1=plt.figure()
plt.hist(loc_prob, bins=[q for q in arange(0,(N+1)*L + bin_width,bin_width)], normed=False, weights=weights, stacked=True)
plt.title('Location probability')
plt.xlabel('pos')
plt.ylabel('#')
plt.grid(True)
plt.xlim(0.0, 1.0)
pp.savefig(fig1)
#plt.clf()

#---------------------------------------------------------------------------------
time_plots=1
time_start=0
time_end=time_start + 100000
markersize=3

scatter_plots=0
scatter_start=555000
scatter_end=565000
#---------------------------------------------------------------------------------

if time_plots==1:
	#---------------------------------------------------------------------------------
	# positions
	y_axis=particle_positions[time_start:time_end]	
	t_axis=xrange(0, len(y_axis))

	fig2=plt.figure(figsize=(20,5))
	plt.title('Trajectory')
	plt.xlabel('t')
	plt.ylabel('pos')
	plt.plot(t_axis,y_axis)#,marker='o',markersize=markersize)
	#plt.axis([time_start,time_end,0,1])
	plt.grid(True)
	pp.savefig(fig2)
	#---------------------------------------------------------------------------------
	# velocities
	y_axis=particle_velocities[time_start:time_end]
	t_axis=xrange(0, len(y_axis))

	fig3=plt.figure(figsize=(20,5))
	plt.title('Velocity (pre)')
	plt.xlabel('t')
	plt.ylabel('v')
	plt.plot(t_axis,y_axis)
	plt.grid(True)
	pp.savefig(fig3)
	#---------------------------------------------------------------------------------
	# lattice velocities pre
	y_axis=lattice_velocities_pre[time_start:time_end]	
	t_axis=xrange(0, len(y_axis))

	fig3=plt.figure(figsize=(20,5))
	plt.title('Current lattice velocity')
	plt.xlabel('t')
	plt.ylabel('xdot')
	plt.plot(t_axis,y_axis)
	plt.grid(True)
	pp.savefig(fig3)
	#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------

if scatter_plots==1:
	#---------------------------------------------------------------------------------
	# recursion plot velocities
	liste1=particle_velocities[:-1]
	liste2=particle_velocities[1:]
	liste1=liste1[scatter_start:scatter_end]
	liste2=liste2[scatter_start:scatter_end]

	fig5=plt.figure()
	plt.title('Recursion plot - velocities (pre)')
	plt.xlabel('v_i')
	plt.ylabel('v_i+1')
	plt.scatter(liste1,liste2,s=0.1)
	plt.grid(True)
	pp.savefig(fig5)
	#---------------------------------------------------------------------------------
	# recursion plot lattice velocities
	liste1=lattice_velocities_pre[:-1]
	liste2=lattice_velocities_pre[1:]
	liste1=liste1[scatter_start:scatter_end]
	liste2=liste2[scatter_start:scatter_end]

	fig6=plt.figure()
	plt.scatter(liste1,liste2,s=0.1)
	plt.title('Recursion plot - velocity of current lattice point (pre)')
	plt.xlabel('xdot_i')
	plt.ylabel('xdot_i+1')
	plt.grid(True)
	pp.savefig(fig6)
	#---------------------------------------------------------------------------------
	# scatter plot mixed velocities
	liste1=particle_velocities
	liste2=lattice_velocities_pre
	liste1=liste1[scatter_start:scatter_end]
	liste2=liste2[scatter_start:scatter_end]

	fig6=plt.figure()
	plt.scatter(liste1,liste2,s=0.1)
	plt.title('Scatter plot - velocity of particle (pre) and current lattice mass')
	plt.xlabel('v_i')
	plt.ylabel('xdot_i')
	plt.grid(True)
	pp.savefig(fig6)
#---------------------------------------------------------------------------------

pp.close()

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
