// particle class file
#ifndef PARTICLE_H
#define PARTICLE_H

#include <math.h>

using namespace std;

const double M=100.0; // particle mass
const double D=0.1; // particle diameter

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
		double* Normed(double*);
		double* Cross(double*, double*, double*);
		double* Cross(double*, int*, double*);
		void Reflect(double*);
		void FoldBack();
		double Touch(double*, int*);
		double* Collide(double, int*, double*);
		double Bounce(double);
		void Evolve(Verlet*);
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

// necessary overload for working with r_touched
double* Particle::Cross(double* vec1, int* vec2, double* res) {
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
	for (int i=0; i<3; i++) {
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

double Particle::Touch(double* r1, int* r2) {
	// see green notebook, 08.12.14
	double* V_normed=new double[3];
	V_normed=this->Normed(V);
	
	double* r2_minus_r1=new double[3]; // include lib or write vector addition
	for (int i=0; i<3; i++) {
		r2_minus_r1[i]=r2[i]-r1[i];
	}
	
	double* r_perp=new double[3];
	for (int i=0; i<3; i++) {
		r_perp[i]=r2[i]-r1[i] - this->Dot(r2_minus_r1, V_normed)*V_normed[i];
	}
	r_perp=this->Normed(r_perp);
	
	// check this
	double lambda=(-1)*this->Dot(r2_minus_r1, r_perp); // perpendicular distance from r2 to bypassing trajectory
	double mu=(-1)*this->Dot(r2_minus_r1, V_normed); // mu gives the position along the velocity line from r1 corresponding to smallest (ie. perpendicular) distance
	
	//cout << lambda << "  " << mu << '\n';
	
	// collision happens if perpendicular distance is smaller than the sum of the two radii
	if (fabs(lambda) < (D/2 + d/2)) return fabs(mu) - sqrt((D+d)*(D+d)/4 - lambda*lambda); // mu should be always positive here
	else return 0;
}

double* Particle::Collide(double m, int* r, double* v) {
	// see green notebook, 08.12.14
	double* p=new double[3];
	double* t=new double[3];
	double* n=new double[3];
	
	for (int i=0; i<3; i++) {
		p[i]=R[i]-r[i];
	}
	p=this->Normed(p);
	
	n=this->Cross(R, r, n);
	t=this->Cross(n, p, t);
	t=this->Normed(t);
	double V_prime=2*(M*this->Dot(V, p) + m*this->Dot(v, p))/(M+m) - this->Dot(V, p); // check this and optimize here
	double v_prime=2*(M*this->Dot(V, p) + m*this->Dot(v, p))/(M+m) - this->Dot(v, p);
	for (int i=0; i<3; i++) {
		VV[i]=V[i]; // do not forget temp cp
		V[i]=this->Dot(VV, t)*t[i] + V_prime*p[i];
		v[i]=this->Dot(v, t)*t[i] + v_prime*p[i];
	}
	
	return v;
}

double Particle::Bounce(double Touched) {
	cout << "touching at: " << Touched << '\n';
	
	// update position
	double* V_normed=new double[3]; // change this later and avoid extra def.
	V_normed=this->Normed(VV);
	for (int i=0; i<3; i++) {
		R[i]=RR[i] + Touched*V_normed[i]; // Is this really pointing into the correct direction?
		RR[i]=R[i]; // Is it okay to update the old positions here?
	}
	
	// get travel time
	double tt=Touched/sqrt(this->Dot(VV, VV));
	return tt;
}

void Particle::Evolve(Verlet* Obj) {
	// determine grid point nearest to particle position
	double* r=new double[3];
	for (int i=0; i<3; i++) {
		r[i]=round(RR[i]);
		cout << RR[i] << '\n';
	}
	
	int X=VV[0]/fabs(VV[0]);
	int Y=VV[1]/fabs(VV[1]);
	int Z; if (fabs(VV[2]) == 0.0) Z=0; else Z=VV[2]/fabs(VV[2]); // variable Z must be declared outside of clause...
	
	cout << '\n';
	
	// distance to border
	int qx; if (X==1) qx=N_X-1 - r[0]; else if (X==-1) qx=r[0]; // perhaps introduce array q...
	int qy; if (Y==1) qy=N_Y-1 - r[1]; else if (Y==-1) qy=r[1];
	int qz; if (Z==1) qz=N_Z-1 - r[2]; else if (Z==-1) qz=r[2]; else if (Z==0) qz=0;
	
	// touching
	double touched=0;
	int* r_touched=new int[3]; // intness causes some trouble...
		// check all grid points in the "direction" of the velocity...
		for (int x=min(0, qx*X); x<max(0+1, qx*X+1); x++) {
			for (int y=min(0, qy*Y); y<max(0+1, qy*Y+1); y++) {
				for (int z=min(0, qz*Z); z<max(0+1, qz*Z+1); z++) {
				
					if ((x==0) && (y==0) && (z==0)) continue; // exclude "self check"
					
					
					r_touched[0]=r[0]+x; r_touched[1]=r[1]+y; r_touched[2]=r[2]+z; // assuming distance between grid points to be unity
					
					touched=this->Touch(RR, r_touched);
					//cout << "touched: " << touched << '\n';	
					if (touched != 0) goto TOUCHDOWN;
					
					//cout << r[0]+x << "  " << r[1]+y << "  " << r[2]+z << "  " << '\n';
					//cout << round(Obj->r0[Obj->Index(r[0]+x, r[1]+y, r[2]+z, 0)]) << '\n'; // for this to work, spatial deplacements of grid points must be small (-> correct rounding)
				}
			}
		}
		TOUCHDOWN: if (touched != 0) {
			// displace particle to next point and evolve grid according to travel time
			Obj->T=this->Bounce(touched);
			Obj->Evolve();
			
			// collision (velocity is being updated inside)
			double* rdot_temp=new double[3]; // at some point I have to get rid of all those 3-arrays
			for (int i=0; i<3; i++) {
				rdot_temp[i]=Obj->rdot[Obj->Index(r_touched[0], r_touched[1], r_touched[2], i)];
			}		
			rdot_temp=this->Collide(m, r_touched, rdot_temp); // Is this the correct and up-to-date velocity?
			
			//change Verlet r1 according to new velocity: r1=dt*rdot + r0 - again: Is this really correct?
			for (int i=0; i<3; i++) {
				int index=Obj->Index(r_touched[0], r_touched[1], r_touched[2], i);
				Obj->rdot[index]=rdot_temp[i];
				Obj->r1[index]=Obj->dt*rdot_temp[i] + Obj->r0[index];
			}
		}
		
		// For tomorrow: Am I now in the position to call Evolve again? Is everything correct and updated?
		
		// How to proceed if grid is left??? reflection etc....
		if (touched==0) cout << "Did not touch!" << '\n';
	
}

#endif
