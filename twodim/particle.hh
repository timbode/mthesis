// particle class file
#ifndef PARTICLE_H
#define PARTICLE_H

#include <math.h>

using namespace std;

const double M=100.0;

class Particle {
	public:
		 Particle(double*, double*); // constructor
		~Particle(); // destructor
		
		double* R; // coordinates
		double* V; // velocities
		double* RR; // previous coordinates
		double* VV; // previous velocities
		
		// methods
		void Reflect(double*);
		void FoldBack();
		double* Collide(double, double*);
};

Particle::Particle(double* R_0, double* V_0) {
	R=new double[3];
	V=new double[3];
	RR=new double[3];
	VV=new double[3];
	for (int i=0; i<3; i++) {
		R[i]=R_0[i];
		V[i]=V_0[i];
		RR[i]=R_0[i];
		VV[i]=V_0[i];
		
	}
}

Particle::~Particle() {

}

void Particle::Reflect(double* n) {
	// normalize normal vector
	double norm_n=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0]=n[0]/norm_n; n[1]=n[1]/norm_n; n[2]=n[2]/norm_n;
	
	// reflect V
	double V_dot_n=V[0]*n[0]+V[1]*n[1]+V[2]*n[2];
	for (int i=0; i<3; i++) {
		V[i]=V[i]-2*V_dot_n*n[i];
	}
}

void Particle::FoldBack() {
	cout << N_X << '\n'; // use dimensions to describe boundary in the box case
	if (R[0] > N_X-1) cout << R[0]-(N_X-1) << '\n';
	if (R[1] > N_Y-1) cout << R[1]-(N_Y-1) << '\n';
	if (R[2] > N_Z-1) cout << R[2]-(N_Z-1) << '\n';
}

double* Particle::Collide(double m, double* v) {
	// see green notebook, 08.12.14
	double V_prime=(M-m)/(M+m);
	double v_prime=2*M/(M+m);
	for (int i=0; i<3; i++) {
		VV[i]=V[i];
		V[i]=V_prime*(VV[i]-v[i]) + v[i];
		v[i]=v_prime*(VV[i]-v[i]) + v[i];
	}
	return v;
}

#endif
