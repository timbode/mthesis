# execution file

num=$(awk '/_DATA_/{ print NR; exit }' constants.hh)
c=1
while read -r line; do
    if [ $c == $num ]; then
        folder=$line
        break
    fi
    c=$((c+1))
done < constants.hh
folder="${folder// /}" # remove spaces
folder="${folder:19:-2}" # extract folder

# make sure that directories exist
mkdir -vp $folder $folder/chunks $folder/hist $folder/init $folder/plots

# number of parallel threads
export OMP_NUM_THREADS=6

# compile
icc -fopenmp -O3 -o main main.cc

#loop particles and repetitions
stats=1
repeat=10
start=0
time {
for ((p=0; p<$stats; ++p)); do
  for ((rep=$start; rep<$repeat; ++rep)); do
      nice -n 19 ./main $repeat $p $rep;
  done
done
echo -e "\n"
}

echo -e "\n========================================================"
echo    "Data evaluation"
echo -e "========================================================\n"

time /home/students/lappet/anaconda/bin/python read_data.py $folder
     /home/students/lappet/anaconda/bin/python make_hist.py $folder
