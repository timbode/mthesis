// Verlet method for a multidimensional grid
#ifndef VERLET_H
#define VERLET_H

#include <vector>
#include <math.h>
#include <omp.h>
#include <algorithm> // max, min

#include "constants.hh"

using namespace std;
using namespace Constants;

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
		void ReadIn();
		int Index(int, int, int, int);
		double NearestNeighbours(int, int, int, int);
		double PotentialEnergy(int, int, int, int, int);
		double Step();
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
	for (int x=0; x<N_[0]; ++x) {
		for (int y=0; y<N_[1]; ++y) {
			for (int z=min(0, N_[2]); z<max(0, N_[2]); ++z) {
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

void Verlet::ReadIn() {
	ifstream state_data("grid_positions.dat");
	int i=0;
	double rr0; double rr1;
	while (state_data >> rr0 >> rr1) {
		r0[i]=rr0; r1[i]=rr1;
		++i;
	}
}

// ------------------------------------------------------------------------------------------------
int Verlet::Index(int X, int Y, int Z, int Alpha) {
	return X*N_[1]*N_[2] + Y*N_[2] + Z + Alpha*N_[0]*N_[1]*N_[2];
}

double Verlet::NearestNeighbours(int X, int Y, int Z, int Alpha) {
	double S=r1[this->Index(X-1, Y, Z, Alpha)] + r1[this->Index(X+1, Y, Z, Alpha)] + r1[this->Index(X, Y-1, Z, Alpha)] + r1[this->Index(X, Y+1, Z, Alpha)];
	if (N_[2]>=3) S+=r1[this->Index(X, Y, Z-1, Alpha)] + r1[this->Index(X, Y, Z+1, Alpha)];
	return S;
}

double Verlet::PotentialEnergy(int Ind, int X, int Y, int Z, int Alpha) {
	int index=Ind;
	double V=0; double term;
	term=r1[index] - r1[this->Index(X-1, Y, Z, Alpha)]; // optimize...
	V+=term*term;
	term=r1[index] - r1[this->Index(X+1, Y, Z, Alpha)];
	V+=term*term;
	term=r1[index] - r1[this->Index(X, Y-1, Z, Alpha)];
	V+=term*term;
	term=r1[index] - r1[this->Index(X, Y+1, Z, Alpha)];
	V+=term*term;
	if (N_[2]>=3) {
		term=r1[index] - r1[this->Index(X, Y, Z-1, Alpha)];
		V+=term*term;
		term=r1[index] - r1[this->Index(X, Y, Z+1, Alpha)];
		V+=term*term;
	}
	return 0.5*k*V;
}

double Verlet::Step() {
		double E; // energy
		#pragma omp parallel for collapse(4) reduction(+:E)
		for (int alpha=0; alpha<3; ++alpha) { 
			for (int x=1; x<(N_[0]-1); ++x) { // Postillon: +++ ++x supposed to be faster than x++ +++
				for (int y=1; y<(N_[1]-1); ++y) {
					for (int z=min(1, N_[2]-1); z<max(1, N_[2]-1); ++z) {
						int index=this->Index(x, y, z, alpha);
						r2[index]=2*r1[index] - r0[index] - ((k/m)*dt*dt)*(2*dim*r1[index] - NearestNeighbours(x, y, z, alpha)); //sqrt ??? ((k/m)*dt*dt) must be << 1
						rdot[index]=(r2[index] - r1[index])/dt; // (r2[index] - r0[index])/(2*dt)
						//cout << r1[index] << "  " << r0[index] << '\n';
						E+=0.5*m*rdot[index]*rdot[index] + PotentialEnergy(index, x, y, z, alpha);
					}
				}
			}
		}
		
				/*for (int alpha=0; alpha<3; ++alpha) { 
			for (int x=0; x<(N_[0]); ++x) { // Postillon: +++ ++x supposed to be faster than x++ +++
				for (int y=0; y<(N_[1]); ++y) {
					for (int z=min(1, N_[2]-1); z<max(1, N_[2]-1); ++z) {
						int index=this->Index(x, y, z, alpha);
						if (x!=0 && y!=0 && x!=N_[0]-1 && y!=N_[1]-1) {
							r2[index]=2*r1[index] - r0[index] - ((k/m)*dt*dt)*(2*dim*r1[index] - NearestNeighbours(x, y, z, alpha)); //sqrt ??? ((k/m)*dt*dt) must be << 1
							rdot[index]=(r2[index] - r1[index])/dt; // (r2[index] - r0[index])/(2*dt)
						}
						
						E+=0.5*m*rdot[index]*rdot[index] + PotentialEnergy(index, x, y, z, alpha);
					}
				}
			}
		}*/

		r0=r1; r1=r2; r2=r0;
		
		return E;
}

void Verlet::Evolve() {
	for (int t=0; t<T/dt; t++) {
		this->Step();
	}
}

// ------------------------------------------------------------------------------------------------
#endif
