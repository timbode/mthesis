// Verlet method for a twodimensional grid
#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>

using namespace std;

const unsigned int dim=2; // dimensions=2,3 
const unsigned int N=3; // Nx...xN grid
const unsigned int Max=pow(N, dim)*dim; // length of arrays
const double L=1.0; // distance between nearest neighbours
const double m=1.0; // mass of grid point
const double kk=1.0; // spring constant (kk in order to avoid ambiguities with indices below)

// ------------------------------------------------------------------------------------------------

class Verlet {
	public:
		Verlet(double, double); // constructor
		Verlet(double, double, double*, double*, double*); // overload constructor
		~Verlet(); // destructor
		
		double T; // time interval
		double dt; // recursion time
		
		double* r0; double* r1; double* r2; // positions: three arrays required because Verlet method depends on the last two sets of positions
		
		double* rdot; // velocities
		
		// methods
		int Index(int, int, int);
		double NearestNeighbours(int, int, int);
		double* Evolve2D();
		int Index(int, int, int, int);
		double NearestNeighbours(int, int, int, int);
		double* Evolve3D();
		
}; // do not forget the funny semicolon down here...

// init
Verlet::Verlet(double T_0, double dt_0) {
	T=T_0; dt=dt_0;
	r0  =new double[Max];
	r1  =new double[Max];
	r2  =new double[Max];
	rdot=new double[Max];
	
	// default initialization: positions must be set
	if (dim==2) {
		for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) {
				  r0[Index(i, j, 0)]=i;   r0[Index(i, j, 1)]=j;
				  r1[Index(i, j, 0)]=i;   r1[Index(i, j, 1)]=j;
				rdot[Index(i, j, 0)]=0; rdot[Index(i, j, 1)]=0;
			}
		}
	}
	if (dim==3) {
		for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) {
				for (int k=0; k<N; k++) {
					  r0[Index(i, j, k, 0)]=i;   r0[Index(i, j, k, 1)]=j;   r0[Index(i, j, k, 2)]=k; // does not he run nicely through the arrays in this way?
					  r1[Index(i, j, k, 0)]=i;   r1[Index(i, j, k, 1)]=j;   r1[Index(i, j, k, 2)]=k;
					rdot[Index(i, j, k, 0)]=0; rdot[Index(i, j, k, 1)]=0; rdot[Index(i, j, k, 2)]=0;
				}
			}
		}
	}
	
}

Verlet::Verlet(double T_0, double dt_0, double* r0_0, double* r1_0, double* rdot_0) {
	T=T_0; dt=dt_0;
	r0  =new double[Max];
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
			if ((j==0) || (j==(N-1))) continue;
			i=(q-alpha-j*N)/(N*N);
			if ((i==0) || (i==(N-1))) continue;
			
			r2[Index(i, j, alpha)]=2*r1[Index(i, j, alpha)] - r0[Index(i, j, alpha)] - (kk*dt*dt/m)*(r1[Index(i, j, alpha)] - NearestNeighbours(i, j, alpha));
			
			// grab velocities
			rdot[Index(i, j, alpha)]=(r2[Index(i, j, alpha)] - r1[Index(i, j, alpha)])/dt;

		}

		/*
		// we leave outer points unchanged in any case (ie. i=1, i<(N-1) etc.) - the square box is then completely defined
		#pragma omp parallel for collapse(4)
		for (int i=1; i<(N-1); i++) {
			for (int j=1; j<(N-1); j++) {
				for (int alpha=0; alpha<dim; alpha++) {
					r2[Index(i, j, alpha)]=2*r1[Index(i, j, alpha)] - r0[Index(i, j, alpha)] - (kk*dt*dt/m)*(r1[Index(i, j, alpha)] - NearestNeighbours(i, j, alpha));
					
					// grab velocities
					rdot[Index(i, j, alpha)]=(r2[Index(i, j, alpha)] - r1[Index(i, j, alpha)])/dt;
				}
			}
		}
		
		*/
		
		
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
			if ((k==0) || (k==(N-1))) continue; // we leave outer points unchanged in any case - the square box is then completely defined
			j=((q-alpha-k*N)/(N*N))%N;
			if ((j==0) || (j==(N-1))) continue;
			i=(q-alpha-k*N-j*N*N)/(N*N*N);
			if ((i==0) || (i==(N-1))) continue;
			
			r2[Index(i, j, k, alpha)]=2*r1[Index(i, j, k, alpha)] - r0[Index(i, j, k, alpha)] - (kk*dt*dt/m)*(r1[Index(i, j, k, alpha)] - NearestNeighbours(i, j, k, alpha));
			
			// grab velocities
			rdot[Index(i, j, k, alpha)]=(r2[Index(i, j, k, alpha)] - r1[Index(i, j, k, alpha)])/dt;
		}
		
		/*
		// we leave outer points unchanged in any case (ie. i=1, i<(N-1) etc.) - the square box is then completely defined
		#pragma omp parallel for collapse(4)
		for (int i=1; i<(N-1); i++) {
			for (int j=1; j<(N-1); j++) {
				for (int k=1; k<(N-1); k++) {
					for (int alpha=0; alpha<dim; alpha++) {
						r2[Index(i, j, k, alpha)]=2*r1[Index(i, j, k, alpha)] - r0[Index(i, j, k, alpha)] - (kk*dt*dt/m)*(r1[Index(i, j, k, alpha)] - NearestNeighbours(i, j, k, alpha));
						// grab velocities
						rdot[Index(i, j, k, alpha)]=(r2[Index(i, j, k, alpha)] - r1[Index(i, j, k, alpha)])/dt;
					}
				}
			}
		}
		
		*/
		
		r0=r1; r1=r2;
		
		return r2;
}
// ------------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------

int main() {



return 0;
}
