# submit to cluster
# /home/students/lappet/cluster
stats=2
for ((p=0; p<$stats; ++p)); do
	mkdir -vp sysfolder_$p
	for s in verlet.hh particle.hh droplet.hh constants.hh main.cc system_data.py read_data.py make_hist.py exec.sh; do
		cp $s sysfolder_$p
	done

	cd sysfolder_$p
	
	# create sub file here
	echo "#$ -N TIMsJob-$p" >> job_$p.sub
	echo "#$ -M lappet@student.ethz.ch" >> job_$p.sub
	echo "#$ -m e" >> job_$p.sub
	echo "#$ -m a" >> job_$p.sub
	echo "#$ -m b" >> job_$p.sub
	echo "#$ -m s" >> job_$p.sub
	# echo "#$ -l vf=100M" >> job_$p.sub
	echo "#$ -R y" >> job_$p.sub
	echo "#$ -pe omp 1" >> job_$p.sub
	echo "nohup ./exec.sh &" >> job_$p.sub
	
	# submit to queue
	#qsub job_$p.sub
	cd ..
done
