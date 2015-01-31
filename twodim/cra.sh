# discard droplets that crashed

stats=3

counter=0
for ((k=0; k<$stats; ++k)); do
	if [ -e sysfolder_$k/data/crashed.dat ]; then
		rm -r sysfolder_$k
		((++counter))
	fi
done
echo -e "$counter directories have been removed!"
