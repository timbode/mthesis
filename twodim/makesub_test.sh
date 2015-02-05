# test
# /home/students/lappet/cluster
stats=1

mkdir -vp data plots

for ((p=0; p<$stats; ++p)); do
	mkdir -vp sysfolder_$p
	for s in verlet.hh particle.hh droplet.hh constants.hh main.cc system_data.py read_data.py make_hist.py exec.sh; do
		cp $s sysfolder_$p
	done

	cd sysfolder_$p

	./exec.sh

	cd ..
done

cp -v sysfolder_0/data/system.dat data/system.dat
