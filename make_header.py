# preparing header for system.cc

# read-in file
lines=open('systems.txt').readlines()

# write to header
counter=0
with open('system_constants.hh', 'w') as g:
	g.write('#ifndef SYSTEM_CONSTANTS_H'+'\n')
	g.write('#define SYSTEM_CONSTANTS_H'+'\n')
	g.write('namespace SystemConstants {'+'\n')
	for line in lines:
		if line.startswith('#'): continue
		if line.startswith('!'): break # separator for different sets of constants in systems.txt
		g.write(line)
		counter+=1
	g.write('}'+'\n'+'#endif'+'\n')

# delete system after usage
open('systems.txt', 'w').writelines(lines[(counter + 2):]) # careful: there must be a blank line after "! ----" in systems.txt for this to work properly
