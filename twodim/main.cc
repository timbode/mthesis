// main file
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>

#include "verlet.hh"
#include "particle.hh"

using namespace std;
using namespace Constants;

int main(int argc, char* argv[]) {

unsigned int repeat=atof(argv[1]);
unsigned int p=atof(argv[2]);
unsigned int rep=atof(argv[3]);

unsigned int steps=1000000;
double dt=1e-5; // should be 1e-5 or 1e-6
double T=steps*dt;

double R_0 [3]={N_[0]/2+0.5, N_[1]/2+0.5, 0.0};
//double V_0 [3]={-0.1, -0.1001, 0};
double V_0 [3]={-1, -1.0001, 0};

/*
Verlet test(T, dt);

test.r1[test.Index(1, 1, 0, 0)]+=0.05;
test.r0[test.Index(1, 1, 0, 0)]+=0.05;
test.r1[test.Index(1, 1, 0, 1)]+=0.05;
test.r0[test.Index(1, 1, 0, 1)]+=0.05;


// open file
ofstream grid_data;
grid_data.open("data/grid.dat");

for (int tt=0; tt<T/dt; tt++) {
	for (int x=0; x<N_[0]; x++) {
		for (int y=0; y<N_[1]; y++) {
			for (int z=0; z<N_[2]; z++) {
				for (int alpha=0; alpha<3; alpha++) {
					grid_data << test.r1[test.Index(x, y, z, alpha)] << ",";
				}
				grid_data << '\t';
			}
		}
	}
	grid_data << '\n';
	test.Evolve();
}
grid_data.close();
*/

unsigned int stats=1;

ofstream system_data;
system_data.open("data/system.dat");

system_data << "# dim: " << dim << '\n';
system_data << "# N_X: " << N_[0] << '\n';
system_data << "# N_Y: " << N_[1] << '\n';
system_data << "# N_Z: " << N_[2] << '\n';
system_data << '\n';
system_data << "# k: " << k << '\n';
system_data << "# m: " << m << '\n';
system_data << "# d: " << d << '\n';
system_data << '\n';
system_data << "# M: " << M << '\n';
system_data << "# D: " << D << '\n';
system_data << '\n';
system_data << "# steps: " << steps << '\n';
system_data << "# repeat: " << repeat << '\n';
system_data << "# dt: " << dt << '\n';
system_data << "# stats: " << stats << '\n';
system_data << '\n';

system_data.close();

//vector< vector<double> > R_0s(stats, vector<double>(3));
//vector< vector<double> > V_0s(stats, vector<double>(3));

		cout << "Starting with repetition " << rep << '\n';
		// create grid instance
		Verlet grid(p, rep, T, dt);
		//grid.r1[grid.Index(1, 1, 0, 0)]+=5;
		//grid.r0[grid.Index(1, 1, 0, 0)]+=5;

		vector<double> data_array((dim+2)*steps, 0.0);

// create particle instance
Particle particle(p, rep, T, dt, R_0, V_0);

// give grid to particle and go
particle.Evolve(&grid, &data_array[0]);

// open file
ofstream particle_data;
ostringstream FileNameStream;
FileNameStream << "data/chunks/particle_" << p << "_chunk_" << rep << ".dat";
string FileName=FileNameStream.str();
particle_data.open(FileName.c_str());

// write array to file
for (int t=0; t<steps; ++t) {
	for (int u=0; u<(dim+2); ++u) {
		particle_data << data_array[(dim+2)*t+u] << '\t';
	}
	particle_data << '\n';
}

particle_data.close();

return 0;
}
