// particle class file
#ifndef PARTICLE_H
#define PARTICLE_H

#include <math.h>

class Particle {
	public:
		 Particle(double*, double*); // constructor
		~Particle(); // destructor
		
		double* r; // coordinates
		double* v; // velocities
		double* rr; // previous coordinates
		double* vv; // previous velocities
		
		// methods
		void Reflect(double*);
		void FoldBack();
};

Particle::Particle(double * r_0, double* v_0) {
	r=new double[3];
	v=new double[3];
	rr=new double[3];
	vv=new double[3];
	for (int i=0; i<3; i++) {
		r[i]=r_0[i];
		v[i]=v_0[i];
		rr[i]=r_0[i];
		vv[i]=v_0[i];
		
	}
}

Particle::~Particle() {

}

void Particle::Reflect(double* n) {
	// normalize normal vector
	double norm_n=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0]=n[0]/norm_n; n[1]=n[1]/norm_n; n[2]=n[2]/norm_n;
	
	// reflect v
	double v_dot_n=v[0]*n[0]+v[1]*n[1]+v[2]*n[2];
	v[0]=v[0]-2*v_dot_n*n[0];
	v[1]=v[1]-2*v_dot_n*n[1];
	v[2]=v[2]-2*v_dot_n*n[2];
}

void Particle::FoldBack() {
	cout << N_X << '\n'; // use dimensions to describe boundary in the box case
	if (r[0] > N_X-1) cout << r[0]-(N_X-1) << '\n';
	if (r[1] > N_Y-1) cout << r[1]-(N_Y-1) << '\n';
	if (r[2] > N_Z-1) cout << r[2]-(N_Z-1) << '\n';
}

#endif
