# execution file

export OMP_NUM_THREADS=6
#icc -openmp -O3 -o system system.cc
g++ -fopenmp -O2 -o system system.cc
nice -n 19 time -f "\n%E real\n%U user\n%S sys\n" ./system
/home/students/lappet/anaconda/bin/python read_system_data.py
