# discard droplets that crashed

stats=500

counter=0
for ((k=0; k<$stats; ++k)); do
	if [ -e sysfolder_$k/data/crashed.dat ]; then
		mv sysfolder_$k crashed_sysfolder_$k
		#rm -r sysfolder_$k
		((++counter))
	fi
done
echo -e "$counter directories have been (re)moved!"
