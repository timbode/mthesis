# execution file

# create location for results
mkdir -vp results

rm -r data
mkdir data

# temporarily copy systems.txt
cp -v systems.txt temp.txt
echo -e "\n"

# specify number of to-be-evaluated systems
for k in 1
	do
		echo -e "\n----- starting with system number $k -----\n"
		../../anaconda/bin/python make_header.py

		export OMP_NUM_THREADS=10
		g++ -fopenmp -O2 -o system system.cc
		nice -n 19 time -f "\n%E real\n%U user\n%S sys\n" ./system
		../../anaconda/bin/python read_system_data.py

		rename data data$k data
		cp -vr data$k results
		rm -vr data$k 
		mkdir -vp data
		echo -e "\n----- system number $k completed -----\n"
	done

# restore systems.txt
cp -v temp.txt systems.txt
rm -v temp.txt

