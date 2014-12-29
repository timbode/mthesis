# execution file

# number of parallel threads
export OMP_NUM_THREADS=6

# compile
icc -fopenmp -O3 -o main main.cc

#loop particles and repetitions
stats=1
repeat=200
start=0
time {
for ((p=0; p<$stats; ++p)); do
  for ((rep=$start; rep<$repeat; ++rep)); do
    nice -n 19 ./main $repeat $p $rep;
  done
done
}

/home/students/lappet/anaconda/bin/python read_data.py
/home/students/lappet/anaconda/bin/python make_hist.py
