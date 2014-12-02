// Verlet method for a twodimensional grid
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

const unsigned int dim=2; // dimensions (can be two or three)
const unsigned int N=3; // Nx...xN grid
const unsigned int Max=pow(N, dim)*dim;
const double L=1.0; // distance between nearest neighbours
const double m=1.0; // mass of grid point
const double K=1.0; // spring constant

class Verlet {
	public:
		Verlet(double, double, double*, double*, double*); // constructor
		~Verlet(); // destructor
		
		double T; // time interval
		double dt; // recursion time
		
		double* r0; double* r1; double* r2; // positions: three arrays required because Verlet method depends on the last two sets of positions
		
		double* rdot; // velocities
		
		// methods
		int Index(int, int, int);
		int Index(int, int, int, int);
		double NearestNeighbours(int, int, int);
		double NearestNeighbours(int, int, int, int);
		double* Evolve2D();
		double* Evolve3D();
		
}; // do not forget the funny semicolon down here...

// init
Verlet::Verlet(double T_0, double dt_0, double* r0_0, double* r1_0, double* rdot_0) {
		T=T_0; dt=dt_0;
		r0  =new double[Max]; // define pow(N, dim)*dim previously
		r1  =new double[Max];
		r2  =new double[Max];
		rdot=new double[Max];
		r0=r0_0; r1=r1_0; rdot=rdot_0;
}

Verlet::~Verlet() {
	
}

// ------------------------------------------------------------------------------------------------
// twodimensional case
int Verlet::Index(int I, int J, int Alpha) {
	return I*N*N + J*N + Alpha;
}

double Verlet::NearestNeighbours(int I, int J, int Alpha) {
	return r1[Index(I-1, J, Alpha)] + r1[Index(I+1, J, Alpha)] + r1[Index(I, J-1, Alpha)] + r1[Index(I, J+1, Alpha)] ;
}

double* Verlet::Evolve2D() {
		int i=0; int j=0;
		int alpha=0;
		
		for (int q=0; q<Max; q++) {
			alpha=q%N;
			j=((q-alpha)/N)%N;
			i=(q-alpha-j*N)/(N*N);
			r2[Index(i, j, alpha)]=2*r1[Index(i, j, alpha)] - r0[Index(i, j, alpha)] - (K*dt/m)*(r1[Index(i, j, alpha)] - NearestNeighbours(i, j, alpha));
			
			// grab velocities
			rdot[Index(i, j, alpha)]=(r2[Index(i, j, alpha)] - r1[Index(i, j, alpha)])/dt;

		}
		
		r0=r1; r1=r2;
		
		return r2;
}
// ------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------
// threedimensional case
int Verlet::Index(int I, int J, int K, int Alpha) {
	return I*N*N*N + J*N*N + K*N + Alpha;
}

double Verlet::NearestNeighbours(int I, int J, int K, int Alpha) {
	return r1[Index(I-1, J, K, Alpha)] + r1[Index(I+1, J, K, Alpha)] + r1[Index(I, J-1, K, Alpha)] + r1[Index(I, J+1, K, Alpha)] + r1[Index(I, J, K-1, Alpha)] + r1[Index(I, J, K+1, Alpha)];
}

double* Verlet::Evolve3D() {
		int i=0; int j=0; int k=0;
		int alpha=0;
		
		for (int q=0; q<Max; q++) {
			alpha=q%N;
			k=((q-alpha)/N)%N; // perhaps optimize here... q-alpha, N*N etc.
			j=((q-alpha-k*N)/(N*N))%N;
			i=(q-alpha-k*N-j*N*N)/(N*N*N);
			r2[Index(i, j, k, alpha)]=2*r1[Index(i, j, k, alpha)] - r0[Index(i, j, k, alpha)] - (K*dt*dt/m)*(r1[Index(i, j, k, alpha)] - NearestNeighbours(i, j, k, alpha));
			
			// grab velocities
			rdot[Index(i, j, k, alpha)]=(r2[Index(i, j, k, alpha)] - r1[Index(i, j, k, alpha)])/dt;
		}
		
		r0=r1; r1=r2;
		
		return r2;
}
// ------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------


int main() {



return 0;
}
