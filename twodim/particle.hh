// particle class file
#ifndef PARTICLE_H
#define PARTICLE_H

#include <math.h>
#include <stdlib.h>

using namespace std;

const double M=100.0; // particle mass
const double D=0.1; // particle diameter

class Particle {
	public:
		 Particle(double, double, double*, double*); // constructor
		~Particle(); // destructor
		
		double T; // time interval
		double dt; // recursion time
		
		double* R; // coordinates
		double* V; // velocities
		
		// methods
		double Dot(double*, double*);
		double* Normed(double*);
		double* Cross(double*, double*, double*);
		void Reflect(double*);
		void FoldBack();
		double Hit(double*, double*);
		double* Collide(double, double*, double*);
		void Evolve(Verlet*);
};

Particle::Particle(double T_0, double dt_0, double* R_0, double* V_0) {
	T=T_0; dt=dt_0;
	R=new double[3];
	V=new double[3];
	for (int i=0; i<3; ++i) {
		R[i]=R_0[i];
		V[i]=V_0[i];
	}
}

Particle::~Particle() {

}

// ------------------------------------------------------------------------------------------------
// Vector analysis

double Particle::Dot(double* vec1, double* vec2) {
	return (vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]);
}

double* Particle::Normed(double* vec) {
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
	n=this->Normed(n);
		
	// reflect V
	double V_dot_n=this->Dot(V, n);
	for (int i=0; i<3; ++i) {
		V[i]=V[i]-2*V_dot_n*n[i];
	}
}

void Particle::FoldBack() {
	cout << N_X << '\n'; // use dimensions to describe boundary in the box case
	if (R[0] > N_X-1) cout << R[0]-(N_X-1) << '\n';
	if (R[1] > N_Y-1) cout << R[1]-(N_Y-1) << '\n';
	if (R[2] > N_Z-1) cout << R[2]-(N_Z-1) << '\n';
	
	// for 2D box
	// specify corner points
	double* r00=new double[3]; r00[0]=0; r00[1]=0; r00[2]=0;
	double* r10=new double[3]; r10[0]=1; r10[1]=0; r10[2]=0;
	double* r01=new double[3]; r01[0]=0; r01[1]=1; r01[2]=0;
	double* r11=new double[3]; r11[0]=1; r11[1]=1; r11[2]=0;
}

// ------------------------------------------------------------------------------------------------
double Particle::Hit(double* r1, double* r2) {
	// see green notebook, 10.12.14
	double* r2_minus_r1=new double[3]; // include lib or write vector addition
	for (int i=0; i<3; ++i) {
		r2_minus_r1[i]=r2[i]-r1[i];
	}
	double l=sqrt(this->Dot(r2_minus_r1, r2_minus_r1));
	return l;
}

double* Particle::Collide(double m, double* r, double* v) {
	// see green notebook, 08.12.14
	double* p=new double[3];
	double* t=new double[3];
	double* n=new double[3];
	
	for (int i=0; i<3; ++i) {
		p[i]=R[i]-r[i];
	}
	p=this->Normed(p);
	
	n=this->Cross(R, r, n);
	t=this->Cross(n, p, t);
	t=this->Normed(t);
	double V_prime=2*(M*this->Dot(V, p) + m*this->Dot(v, p))/(M+m) - this->Dot(V, p); // check this and optimize here
	double v_prime=2*(M*this->Dot(V, p) + m*this->Dot(v, p))/(M+m) - this->Dot(v, p);
	
	double* V_temp=new double[3]; // do not forget temp cp here!
	for (int i=0; i<3; ++i) V_temp[i]=V[i];
	double* v_temp=new double[3]; // do not forget temp cp here!
	for (int i=0; i<3; ++i) v_temp[i]=v[i];
	
	for (int i=0; i<3; ++i) {
		V[i]=this->Dot(V_temp, t)*t[i] + V_prime*p[i];
		v[i]=this->Dot(v_temp, t)*t[i] + v_prime*p[i];
	}
	
	return v;
}

void Particle::Evolve(Verlet* Obj) {
	for (int t=0; t<T/dt; t++) {
		// underlying assumption: displacement of grid points is small enough such that only collisions with the nearest grid point actually occur
		// determine grid point nearest to particle position
		double* r=new double[3];
		for (int i=0; i<3; ++i) r[i]=round(R[i]);

		// determine actual position and velocity
		double* r_nearest=new double[3];
		double* rdot_nearest=new double[3];
		for (int i=0; i<3; ++i) {
			int index=Obj->Index(r[0], r[1], r[2], i);
			r_nearest[i]=Obj->r1[index];
			rdot_nearest[i]=Obj->rdot[index];
			//cout << Obj->r1[index] << ", ";
		}
		cout << '\n';
		
		// check distance to actual position
		if (this->Hit(R, r_nearest) <= (D/2 + d/2)) {
			cout << "Hit!" << '\n';
			// make collision
			rdot_nearest=this->Collide(m, r_nearest, rdot_nearest);
		
			//update r1 according to new velocity: r1=dt*rdot + r0
			for (int i=0; i<3; ++i) {
				int index=Obj->Index(r[0], r[1], r[2], i);
				Obj->rdot[index]=rdot_nearest[i];
				Obj->r1[index]=dt*rdot_nearest[i] + Obj->r0[index];
				//cout  << Obj->r1[index] << ", ";
			}
			cout << '\n';
		}
	
		// evolve grid and particle
		cout << "Particle: ";
		for (int i=0; i<3; ++i) {
			R[i]=R[i] + dt*V[i];
			
			//cout << R[i] << ", ";
		}
		cout << '\n';
		Obj->Step();
	}
}
// ------------------------------------------------------------------------------------------------
#endif
