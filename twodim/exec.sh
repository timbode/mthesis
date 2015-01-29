# execution file

folder="data"

# make sure that directories exist
mkdir -vp $folder $folder/chunks $folder/hist $folder/init $folder/plots

# number of parallel threads
export OMP_NUM_THREADS=6

# compile
icc -openmp -O3 -o main main.cc
#g++ -Wall -fopenmp -O3 -o main main.cc

#loop repetitions
repeat=1
start=0
time {
for ((rep=$start; rep<$repeat; ++rep)); do
          ./main $repeat 0 $rep;
done
echo -e "\n"
}

echo -e "\n========================================================"
echo    "Data evaluation"
echo -e "========================================================\n"

time /home/students/lappet/anaconda/bin/python read_data.py $folder
     /home/students/lappet/anaconda/bin/python make_hist.py $folder
