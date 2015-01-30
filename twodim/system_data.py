# read system.dat
import sys

def read():
    global Dict; Dict={}
    Dict["_DATA_"]=sys.argv[1]

    constants=[]
    with open(Dict["_DATA_"]+'/system.dat') as g:
    	for line in g:
    		if line.startswith('#'):
    			strang=line.strip().split()
    			constants.append(strang)
    	g.close()

    Dict["dim"]=int(constants[0+1][2])
    Dict["N_X"]=int(constants[1+1][2])
    Dict["N_Y"]=int(constants[2+1][2])
    Dict["k"]=float(constants[4+1][2])
    Dict["m"]=float(constants[5+1][2])
    Dict["d"]=float(constants[6+1][2])
    Dict["L"]=float(constants[7+1][2])
    Dict["M"]=float(constants[7+1+1][2])
    Dict["D"]=float(constants[8+1+1][2])
    Dict["f"]=float(constants[9+1+1][2])

    Dict["steps"]=int(constants[9+1+1+1][2])
    Dict["repeat"]=int(constants[10+1+1+1][2])
    Dict["dt"]=float(constants[11+1+1+1][2])
    Dict["T"]=Dict["steps"]*Dict["repeat"]*Dict["dt"]

    Dict["stats"]=int(constants[12+1+1+1][2])
    Dict["system_type"]=str(constants[13+1+1+1][2])

    global bins_x, bins_y
    bins_y=30; bins_x=bins_y*(Dict["N_X"]-1)/(Dict["N_Y"]-1);
