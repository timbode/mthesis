# make historgram
from numpy import *
from math import *
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

matplotlib.rcParams["agg.path.chunksize"]=20000

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

constants=[]
with open('data/system.dat') as f:
	for line in f:
		if line.startswith('#'):
			strang=line.strip().split()
			constants.append(strang)
	f.close()

print "================================================="
print constants
print "================================================="

dim=int(constants[0][2])
N_X=int(constants[1][2])
N_Y=int(constants[2][2])
k=float(constants[4][2])
m=float(constants[5][2])
d=float(constants[6][2])
M=float(constants[7][2])
D=float(constants[8][2])
steps=int(constants[9][2])
repeat=int(constants[10][2])
dt=float(constants[11][2])
stats=int(constants[12][2])
T=steps*repeat*dt

bins_x=50; bins_y=50

# 2D histogram
all_counts=zeros((bins_x-1, bins_y-1))
for p in xrange(0, stats):
	for rep in xrange(0, repeat):
		with open('data/hist/hist_counts_'+str(p)+'_chunk_'+str(rep)+'.dat') as f:
			counts=[]
			for line in f:
				strang=line.strip().split()
				strang=[float(q) for q in strang]
				counts.append(strang)
			counts=transpose(array(counts))
			all_counts+=counts

# binning
xbins=linspace(0, N_X, bins_x)
ybins=linspace(0, N_Y, bins_y)
xbins, ybins=meshgrid(xbins, ybins)

fig=plt.figure()
size=0.7
if N_X > N_Y:
	size_x=size
	size_y=size*(N_Y-1)/(N_X-1)
else:
	size_x=size*(N_X-1)/(N_Y-1)
	size_y=size
ax=fig.add_axes((0.15, 0.2, size_x, size_y))
plt.title('Location probability')
mesh=ax.pcolormesh(xbins, ybins, all_counts)
fig.colorbar(mesh)
font_size=8
x_text=0.1
fig.text(x_text, 0.04, 'Verlet: '+'dt='+str(dt)+', '+'T='+str(T), fontsize=font_size)
fig.text(x_text, 0.07, 'Grid: '+'m='+str(m)+', '+'d='+str(d)+', '+'k='+str(k), fontsize=font_size)
fig.text(x_text, 0.1, 'Particle: '+'M='+str(M)+', '+'D='+str(D), fontsize=font_size)
plt.xlim(0, N_X-1)
plt.ylim(0, N_Y-1)
plt.xlabel('x')
plt.ylabel('y')
#fig.tight_layout()
fig.savefig('data/plots/histogram2D.png')
