# submit group of jobs to cluster
# /home/students/lappet/cluster
stats=1
for ((p=0; p<$stats; ++p)); do
    # create sub file here
    echo "#$ -N tjob-$p" >> job_$p.sub
    echo "#$ -M lappet@student.ethz.ch" >> job_$p.sub
    echo "#$ -m e" >> job_$p.sub
    echo "#$ -m a" >> job_$p.sub
    #echo "#$ -m b" >> job_$p.sub
    echo "#$ -m s" >> job_$p.sub
    echo "#$ -l vf=300M" >> job_$p.sub # include this for the IBM cluster
    echo "#$ -R y" >> job_$p.sub
    #echo "#$ -pe omp 1" >> job_$p.sub
    echo "./groupjobs.sh $p" >> job_$p.sub

    # submit to queue
    qsubn job_$p.sub # qsubn for the IBM cluster
done
