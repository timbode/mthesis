// main file
#include <iostream>
#include <fstream>

#include "verlet.hh"
#include "particle.hh"

using namespace std;

int main() {

int time_steps=1000;
double T=2.0;
double dt=0.01;

Verlet grid(T, dt);

grid.r1[grid.Index(1, 1, 0, 0)]+=0.0;
grid.r0[grid.Index(1, 1, 0, 0)]+=0.0;
grid.r1[grid.Index(1, 1, 0, 1)]+=0.0;
grid.r0[grid.Index(1, 1, 0, 1)]+=0.0;

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

Verlet grid2(T, dt);

double R_0 [3]={1.5,1.5,0.0};
double V_0 [3]={3.00001,3,0};
double n [3]={0,1,0};
Particle particle(T, dt, R_0, V_0);
//particle.Reflect(n);
//cout << particle.V[0] << "   " << particle.V[1] << "   " << particle.V[2] << "   " << '\n';
//particle.FoldBack();
particle.Evolve(&grid2);

/*
cout << '\n';
double* r; r=new double[3];
r[0]=0; r[1]=3.0001; r[2]=3.0;
double* v; v=new double[3];
v[0]=0; v[1]=0.0; v[2]=1.0;
//v=particle.Collide(m, r, v);
//cout << particle.V[0] << "   " << particle.V[1] << "   " << particle.V[2] << "   " << '\n';
//cout << v[0] << "   " << v[1] << "   " << v[2] << "   " << '\n';
cout << particle.Touch(R_0, r) << '\n';
*/
return 0;
}
