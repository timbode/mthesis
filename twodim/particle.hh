// particle class file
#ifndef PARTICLE_H
#define PARTICLE_H

#include <math.h>

using namespace std;

const double M=100.0;
const double D=1.0; // diameter of particle

class Particle {
	public:
		 Particle(double*, double*); // constructor
		~Particle(); // destructor
		
		double* R; // coordinates
		double* V; // velocities
		double* RR; // previous coordinates
		double* VV; // previous velocities
		
		// methods
		double Dot(double*, double*);
		double* Norm(double*);
		double* Cross(double*, double*, double*);
		void Reflect(double*);
		void FoldBack();
		double* Collide(double, double*, double*);
		double Touch(double*, double*);
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

// ------------------------------------------------------------------------------------------------
// Vector analysis

double Particle::Dot(double* vec1, double* vec2) {
	return (vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]);
}

double* Particle::Norm(double* vec) {
	double norm_vec=sqrt(this->Dot(vec, vec));
	vec[0]=vec[0]/norm_vec; vec[1]=vec[1]/norm_vec; vec[2]=vec[2]/norm_vec;
	return vec;
}

double* Particle::Cross(double* vec1, double* vec2, double* res) {
	res[0]=vec1[1]*vec2[2] - vec1[2]*vec2[1];
	res[1]=vec1[2]*vec2[0] - vec1[0]*vec2[2];
	res[2]=vec1[0]*vec2[1] - vec1[1]*vec2[0];
	return res;
}
// ------------------------------------------------------------------------------------------------

void Particle::Reflect(double* n) {
	// normalize normal vector
	n=this->Norm(n);
		
	// reflect V
	double V_dot_n=this->Dot(V, n);
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

double* Particle::Collide(double m, double* r, double* v) {
	// see green notebook, 08.12.14
	double* p=new double[3];
	double* t=new double[3];
	double* n=new double[3];
	
	for (int i=0; i<3; i++) {
		p[i]=R[i]-r[i];
	}
	p=this->Norm(p);
	
	n=this->Cross(R, r, n);
	t=this->Cross(n, p, t);
	t=this->Norm(t);
	double V_prime=2*(M*this->Dot(V, p) + m*this->Dot(v, p))/(M+m) - this->Dot(V, p); // check this
	double v_prime=2*(M*this->Dot(V, p) + m*this->Dot(v, p))/(M+m) - this->Dot(v, p);
	for (int i=0; i<3; i++) {
		VV[i]=V[i]; // do not forget temp cp
		V[i]=this->Dot(VV, t)*t[i] + V_prime*p[i];
		v[i]=this->Dot(v, t)*t[i] + v_prime*p[i];
	}
	
	return v;
}

double Particle::Touch(double* r1, double* r2) {
	// see green notebook, 08.12.14
	double* V_normed=new double[3];
	V_normed=this->Norm(V);
	
	double* r2_minus_r1=new double[3]; // include lib or write vector addition
	for (int i=0; i<3; i++) {
		r2_minus_r1[i]=r2[i]-r1[i];
	}
	
	double* r_perp=new double[3];
	for (int i=0; i<3; i++) {
		r_perp[i]=r2[i]-r1[i] - this->Dot(r2_minus_r1, V_normed)*V_normed[i];
	}
	r_perp=this->Norm(r_perp);
	
	// check this
	double lambda=(-1)*this->Dot(r2_minus_r1, r_perp); // perpendicular distance from r2 to bypassing trajectory
	double mu=(-1)*this->Dot(r2_minus_r1, V_normed); // mu gives the position along the velocity line from r1 corresponding to smallest (ie. perpendicular) distance
	
	// collision happens if perpendicular distance is smaller than the sum of the two radii
	if (lambda < (D/2 + d/2)) return mu - sqrt((D+d)*(D+d)/4 - lambda*lambda); // mu should turn out to be always positive here
	else return 0;
}

#endif
