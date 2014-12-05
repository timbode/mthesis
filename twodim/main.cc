// main file
#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>
#include <fstream>
#include <algorithm>

#include "verlet.hh"

using namespace std;

int main() {

int time_steps=1000;
double T=0.01;
double dt=0.0001;

Verlet grid(T, dt);

grid.r1[grid.Index(1, 1, 0, 0)]+=0.01;
grid.r0[grid.Index(1, 1, 0, 0)]+=0.01;
grid.r1[grid.Index(1, 1, 0, 1)]+=0.01;
grid.r0[grid.Index(1, 1, 0, 1)]+=0.01;

// open file
ofstream grid_data;
grid_data.open("data/grid.dat");

for (int tt=0; tt<time_steps; tt++) {
	for (int x=0; x<N_X; x++) {
		for (int y=0; y<N_Y; y++) {
			for (int z=0; z<N_Z; z++) {
				for (int alpha=0; alpha<3; alpha++) {
					grid_data << grid.r0[grid.Index(x, y, z, alpha)] << ",";					
				}
				grid_data << '\t';
			}
		}
	}
	grid_data << '\n';
	grid.Evolve();
}
grid_data.close();

return 0;
}
