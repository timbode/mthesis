// Verlet method for a multidimensional grid
#ifndef VERLET_H
#define VERLET_H

#include <vector>
#include <math.h>
#include <omp.h>
#include <algorithm> // max, min

using namespace std;

const unsigned int dim=2;
// careful: N_X, N_Y must be >=3; N_Z can be either =1 or >=3 (there have to be inner points for the nearest-neighbour concept to work properly)
const int N_X=3; const int N_Y=3; const int N_Z=1;
const unsigned int Max=3*N_X*N_Y*N_Z; // length of arrays
const double m=100.0; // mass of grid point
const double k=1e3; // spring constant
const double d=0.1; // this is now rather the diameter of the oscillation region than the grid point diameter

// ------------------------------------------------------------------------------------------------

class Verlet {
	public:
		Verlet(double, double); // constructor
		~Verlet(); // destructor
		
		double T; // time interval
		double dt; // recursion time
		
		double* r0; double* r1; double* r2; // positions: three arrays required because Verlet method depends on the last two sets of positions
		
		double* rdot; // velocities
		
		// methods
		int Index(int, int, int, int);
		double NearestNeighbours(int, int, int, int);
		void Step();
		void Evolve();
		
}; // do not forget the funny semicolon down here...

// init
Verlet::Verlet(double T_0, double dt_0) {
	T=T_0; dt=dt_0;
	r0  =new double[Max];
	r1  =new double[Max];
	r2  =new double[Max];
	rdot=new double[Max];
	
	// default initialization: positions must be set
	for (int x=0; x<N_X; ++x) {
		for (int y=0; y<N_Y; ++y) {
			for (int z=min(0, N_Z); z<max(0, N_Z); ++z) {
				  r0[this->Index(x, y, z, 0)]=x;   r0[this->Index(x, y, z, 1)]=y;   r0[this->Index(x, y, z, 2)]=z;
				  r1[this->Index(x, y, z, 0)]=x;   r1[this->Index(x, y, z, 1)]=y;   r1[this->Index(x, y, z, 2)]=z;
				  r2[this->Index(x, y, z, 0)]=x;   r2[this->Index(x, y, z, 1)]=y;   r2[this->Index(x, y, z, 2)]=z;
				rdot[this->Index(x, y, z, 0)]=0; rdot[this->Index(x, y, z, 1)]=0; rdot[this->Index(x, y, z, 2)]=0;
			}
		}
	}
	
}

Verlet::~Verlet() {
	
}

// ------------------------------------------------------------------------------------------------
int Verlet::Index(int X, int Y, int Z, int Alpha) {
	return X*N_Y*N_Z + Y*N_Z + Z + Alpha*N_X*N_Y*N_Z;
}

double Verlet::NearestNeighbours(int X, int Y, int Z, int Alpha) {
	double S=r1[this->Index(X-1, Y, Z, Alpha)] + r1[this->Index(X+1, Y, Z, Alpha)] + r1[this->Index(X, Y-1, Z, Alpha)] + r1[this->Index(X, Y+1, Z, Alpha)];
	if (N_Z>=3) S+=r1[this->Index(X, Y, Z-1, Alpha)] + r1[this->Index(X, Y, Z+1, Alpha)];
	return S;
}

void Verlet::Step() {
		//#pragma omp parallel for collapse(4)
		for (int alpha=0; alpha<3; ++alpha) { 
			for (int x=1; x<(N_X-1); ++x) { // ++x supposed to be faster than x++
				for (int y=1; y<(N_Y-1); ++y) {
					for (int z=min(1, N_Z-1); z<max(1, N_Z-1); ++z) {
						int index=this->Index(x, y, z, alpha);
						r2[index]=2*r1[index] - r0[index] - (k*dt*dt/m)*(2*dim*r1[index] - NearestNeighbours(x, y, z, alpha));
						
						rdot[index]=(r2[index] - r0[index])/(2*dt);
					}
				}
			}
		}

		r0=r1; r1=r2; r2=r0;
}

void Verlet::Evolve() {
	for (int t=0; t<T/dt; t++) { // t=0???
		this->Step();
	}
}

// ------------------------------------------------------------------------------------------------
#endif
