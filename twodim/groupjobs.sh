# group jobs

group_size=10

for ((p=$1*$group_size; p<$1*$group_size+$group_size; ++p)); do
    mkdir -vp sysfolder_$p
    for s in verlet.hh particle.hh droplet.hh constants.hh main.cc system_data.py read_data.py make_hist.py exec.sh; do
        cp $s sysfolder_$p
    done

    cd sysfolder_$p

    ./exec.sh

    cd ..
done
