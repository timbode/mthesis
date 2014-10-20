from numpy import *
from math import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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

burn_in=70000

loc_prob1=[]
weights1=[]
loc_prob2=[]
weights2=[]
loc_prob3=[]
weights3=[]
with open('particle_data.txt') as f:
	counter=0
	initial=0
	for k, line in enumerate(f):
		if line.startswith("#"):
			strang=line.strip().split()
			constants.append(strang)
		else:	
			if k < (counter*int(constants[5][2]) + burn_in + 9):
				continue
				
			if ((k-9) % int(constants[5][2])) == 0:
				counter+=1
				continue
			print k
				
			string=line.strip().split()
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

pp=PdfPages('output.pdf')		
plt.figure()
plt.hist([loc_prob1,loc_prob2,loc_prob3], bins=[q for q in arange(0,(N+1)*L + bin_width,bin_width)], normed=False, weights=[weights1,weights2,weights3], stacked=True)
#plt.ylim(0.0, 0.7e+5)
plt.savefig(pp,format='pdf')
plt.clf()
pp.close()
