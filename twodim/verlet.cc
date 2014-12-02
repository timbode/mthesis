// Verlet method for a twodimensional grid
#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>
#include <fstream>

using namespace std;

const unsigned int dim=2; // dimensions=2,3 
const unsigned int N=4; // (N-1)x...x(N-1) grid for (N-1) > dim - outer points are held fixed, so the effective grid size is reduced
const unsigned int Max=pow(N, dim)*dim; // length of arrays
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
		void Step2D();
		void Evolve2D();
		int Index(int, int, int, int);
		double NearestNeighbours(int, int, int, int);
		void Step3D();
		void Evolve3D();
		
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
				  r2[Index(i, j, 0)]=i;   r2[Index(i, j, 1)]=j;
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
					  r2[Index(i, j, k, 0)]=i;   r2[Index(i, j, k, 1)]=j;   r2[Index(i, j, k, 2)]=k;
					rdot[Index(i, j, k, 0)]=0; rdot[Index(i, j, k, 1)]=0; rdot[Index(i, j, k, 2)]=0;
				}
			}
		}
	}
	
}

Verlet::Verlet(double T_0, double dt_0, double* r0_0, double* r1_0, double* rdot_0) {
	T=T_0; dt=dt_0;
	//r0  =new double[Max];
	//r1  =new double[Max];
	//r2  =new double[Max];
	//rdot=new double[Max];
	r0=r0_0; r1=r1_0; r2=r1_0; rdot=rdot_0;
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

void Verlet::Step2D() {
		// we leave outer points unchanged in any case (therefore i=1, i<(N-1) etc.) - the square box is then completely defined
		//#pragma omp parallel for collapse(3)
		for (int i=1; i<(N-1); i++) {
			for (int j=1; j<(N-1); j++) {
				for (int alpha=0; alpha<dim; alpha++) {
					int index=Index(i, j, alpha);
					r2[index]=2*r1[index] - r0[index] - (kk*dt*dt/m)*(r1[index] - NearestNeighbours(i, j, alpha));
					
					// grab velocities
					rdot[index]=(r2[index] - r1[index])/dt;
				}
			}
		}
		
		r0=r1; r1=r2;
}

void Verlet::Evolve2D() {
	for (int t=0; t<T/dt; t++) {
		this->Step2D();
	}
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

void Verlet::Step3D() {
		#pragma omp parallel for collapse(4)
		for (int i=1; i<(N-1); i++) {
			for (int j=1; j<(N-1); j++) {
				for (int k=1; k<(N-1); k++) {
					for (int alpha=0; alpha<dim; alpha++) {
						int index=Index(i, j, alpha);
						r2[index]=2*r1[index] - r0[index] - (kk*dt*dt/m)*(r1[index] - NearestNeighbours(i, j, k, alpha));
						
						rdot[index]=(r2[index] - r1[index])/dt;
					}
				}
			}
		}
		
		r0=r1; r1=r2;
}

void Verlet::Evolve3D() {
	for (int t=0; t<T/dt; t++) {
		this->Step3D();
	}
}
// ------------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------

int main() {

int time_steps=10;
double T=0.05;

Verlet grid(0.05, 0.0001);

for (int i=0; i<N; i++) {
	for (int j=0; j<N; j++) {
		for (int alpha=0; alpha<dim; alpha++) {
			cout << grid.r2[grid.Index(i, j, alpha)] << ", ";					
		}
	}
}
cout << '\n';
//grid.rdot[grid.Index(1, 1, 0)]=1.0;
//grid=Verlet(0.05, 0.0001, grid.r1, grid.r2, grid.rdot);

// open file
ofstream grid_data;
grid_data.open("data/grid_data.txt");

for (int tt=0; tt<time_steps; tt++) {
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			for (int alpha=0; alpha<dim; alpha++) {
				//grid_data << grid.r2[grid.Index(i, j, alpha)] << ",";					
			}
			grid_data << '\t';
		}
	}
	grid_data << '\n';
	//grid.Evolve2D();
	//grid=Verlet(0.05, 0.0001, grid.r1, grid.r2, grid.rdot);	
}
grid_data.close();
/*
Verlet grid(0.05, 0.0001);
grid.Evolve2D();
grid=Verlet(0.05, 0.0001, grid.r1, grid.r2, grid.rdot);
grid.Evolve2D();
for (int i=0; i<N; i++) {
	for (int j=0; j<N; j++) {
		for (int alpha=0; alpha<dim; alpha++) {
			cout << grid.r0[grid.Index(i, j, alpha)] << ", ";					
		}
	}
}
cout << '\n';

for (int i=0; i<N; i++) {
	for (int j=0; j<N; j++) {
		for (int alpha=0; alpha<dim; alpha++) {
			cout << grid.r2[grid.Index(i, j, alpha)] << ", ";					
		}
	}
}
cout << '\n';
*/
return 0;
}
