// droplet class file
#ifndef DROPLET_H
#define DROPLET_H

#include <math.h>
#include <stdlib.h>
#include <algorithm> // max, min

#include "constants.hh"

using namespace std;
using namespace Constants;

// ------------------------------------------------------------------------------------------------

class Droplet {
	public:
		 Droplet(int, int, unsigned int, double, double*, double*); // constructor
		~Droplet(); // destructor

		int p; // droplet number
		int rep; // repetition number

		unsigned int Steps; // number of steps
		double dt; // recursion time
		unsigned int T_col; // collision time

		double* R; // coordinates
		double* V; // velocities

		// methods
		double Dot(double*, double*);
		double* Normed(double*);
		double* Cross(double*, double*, double*);
		void Reflect(double*);
		bool Hit(double*, double*);
		bool Hit(unsigned int);
		double* Collide(double, double*, double*);
		void Evolve(Verlet*, double*);
};

Droplet::Droplet(int P, int Rep, unsigned int Steps_0, double dt_0, double* R_0, double* V_0) {
	p=P; rep=Rep;
	Steps=Steps_0; dt=dt_0;
	T_col=(1/f)*(1/dt);
	R=new double[3];
	V=new double[3];

	if (rep==0) {
		for (int i=0; i<3; ++i) {
			R[i]=R_0[i];
			V[i]=V_0[i];
		}
	}
	else {
		ifstream state_data;
		ostringstream FileNameStream;
		FileNameStream << _DATA_ << "/init/" << system_type << "_" << p << "_init_chunk_" << rep << ".dat";
		string FileName=FileNameStream.str();
		state_data.open(FileName.c_str());
		int i=0;
		double RR; double VV;
		while (state_data >> RR >> VV) {
			R[i]=RR; V[i]=VV;
			//cout << setprecision(15) << "construction: " << V[i] << "   " << VV << "\n";
			++i;
		}

		state_data.close();
	}
}

Droplet::~Droplet() {
	// save current state of the system
	ofstream state_data;
	ostringstream FileNameStream;
	FileNameStream << _DATA_ << "/init/" << system_type << "_" << p << "_init_chunk_" << rep + 1 << ".dat";
	string FileName=FileNameStream.str();
	state_data.open(FileName.c_str());
	state_data.precision(15);
	// write droplet state to file
	for (int i=0; i<3; ++i) {
		state_data << R[i] << '\t';
		state_data << V[i] << '\t';
		state_data << '\n';
	}

	state_data.close();
}

// ------------------------------------------------------------------------------------------------
// Vector analysis

double Droplet::Dot(double* vec1, double* vec2) {
	return (vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]);
}

double* Droplet::Normed(double* vec) {
	double norm_vec=sqrt(this->Dot(vec, vec));
	for (int i=0; i<3; ++i) vec[i]=vec[i]/norm_vec;
	return vec;
}

double* Droplet::Cross(double* vec1, double* vec2, double* res) {
	res[0]=vec1[1]*vec2[2] - vec1[2]*vec2[1];
	res[1]=vec1[2]*vec2[0] - vec1[0]*vec2[2];
	res[2]=vec1[0]*vec2[1] - vec1[1]*vec2[0];
	return res;
}
// ------------------------------------------------------------------------------------------------

void Droplet::Reflect(double* n) {
	// normalize normal vector
	n=this->Normed(n);

	// reflect V
	double V_dot_n=this->Dot(V, n);
	for (int i=0; i<3; ++i) {
		V[i]=V[i]-2*V_dot_n*n[i];
	}
}

// ------------------------------------------------------------------------------------------------
bool Droplet::Hit(double* r1, double* r2) {
	// see green notebook, 10.12.14
	double* r2_minus_r1=new double[3]; // include lib or write vector addition
	for (int i=0; i<3; ++i) {
		r2_minus_r1[i]=r2[i]-r1[i];
	}
	double l=sqrt(this->Dot(r2_minus_r1, r2_minus_r1));
	return l <= (D/2 + d/2);
}

bool Droplet::Hit(unsigned int tt) {
	return (tt % T_col) == 0;
}

// this function is confirmed to be working properly; set 3x3 grid for inspection
double* Droplet::Collide(double m, double* r, double* v) {
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

	double V_Dot_p=this->Dot(V, p);
	double v_Dot_p=this->Dot(v, p);
	double c_factor=2*(M*V_Dot_p + m*v_Dot_p)/(M+m);

	double V_Dot_t=this->Dot(V, t);
	double v_Dot_t=this->Dot(v, t);

	for (int i=0; i<3; ++i) {
		V[i]=V_Dot_t*t[i] + (c_factor - V_Dot_p)*p[i];
		v[i]=v_Dot_t*t[i] + (c_factor - v_Dot_p)*p[i];
	}

	return v;
}

// underlying assumption: displacement of grid points is small enough such that only collisions with the nearest grid point actually occur
void Droplet::Evolve(Verlet* Obj, double* datarr) {
	for (unsigned int t=0; t<Steps; t++) {
		double E; // droplet energy
		double E_grid; // grid energy

		// determine grid point nearest to droplet position
		double* r=new double[3];
		for (int i=0; i<3; ++i) {
			r[i]=round(R[i]/L); // NOTE: factor of 1/L
			}

		// exlcude collisions with outer grid points and make droplet stay in the box
		double* n=new double[3];
		if ((r[0]==N_[0]-1 || r[0]==0) || (r[1]==N_[1]-1 || r[1]==0) || ((N_[2]!=1) && (r[2]==N_[2]-1 || r[2]==0))) {
			for (int i=0; i<3; ++i) n[i]=min(L*(N_[i]-1) - R[i], R[i]); // n[2] is always zero... NOTE: factor of L
			if ((N_[2]==1) && (n[2]==0)) n[2]=max(n[0], n[1]) + 1; // just make n[2] the biggest
			double s=min(min(n[0], n[1]), n[2]);
			for (int i=0; i<3; ++i) {
				double temp=n[i];
				n[i]=0;
				if (temp == s) {
					n[i]=(1 - (N_[i]-1 - r[i])*2/(N_[i]-1))*fabs(temp); // orientation outwards
				}
			}

			// reflect
			if ((s <= 0) && (this->Dot(n, V) > 0)) {
				this->Reflect(n); // second condition is to avoid that droplet gets stuck in the corner
			}

			// evolve droplet
			for (int i=0; i<3; ++i) {
				R[i]=R[i] + dt*V[i];
				E+=0.5*M*V[i]*V[i]; // energy

				if (i<dim) {
					*datarr=R[i];
					++datarr;
				}
				//cout << V[i] << ", ";
				//cout << R[i] << ", ";
			}
			*datarr=E; ++datarr;
			//cout << '\n';

			// evolve grid
			E_grid=Obj->Step();
			*datarr=E_grid; ++datarr;

			continue;
		}

		// ------------------------------------------------------------------------------------------------

		// determine actual position and velocity
		double* r_nearest=new double[3];
		double* rdot_nearest=new double[3];
		for (int i=0; i<3; ++i) {
			//cout << "trouble: " << r[i] << ", ";
			int index=Obj->Index(r[0], r[1], r[2], i);
			//cout << "index: " << index << "   " <<  Obj->r1[index] << '\n';
			r_nearest[i]=Obj->r1[index];
			rdot_nearest[i]=Obj->rdot[index];
			//cout << "r_nearest: " << r_nearest[i] << ", ";
			//cout << R[i] << ", ";
		}
		//cout << '\n';

		// check if it is time
		if (this->Hit(t)) { // if (this->Hit(R, r_nearest)) {
			//cout << "Hit!" << '\n';
			// make collision
			rdot_nearest=this->Collide(m, r_nearest, rdot_nearest);

			// update r1 according to new velocity: r1=dt*rdot + r0
			for (int i=0; i<3; ++i) {
				int index=Obj->Index(r[0], r[1], r[2], i);
				Obj->rdot[index]=rdot_nearest[i];
				Obj->r1[index]=dt*rdot_nearest[i] + Obj->r0[index];
				//cout  << Obj->rdot[index] << ", ";
			}
			//cout << '\n';
		}

		// evolve droplet
		//cout << "Droplet: ";
		for (int i=0; i<3; ++i) {
			R[i]=R[i] + dt*V[i];
			E+=0.5*M*V[i]*V[i]; // energy

			if (i<dim) {
				*datarr=R[i];
				++datarr;
			}
			//cout << V[i] << ", ";
		}
		*datarr=E; ++datarr;
		//cout << '\n';

		// evolve grid
		E_grid=Obj->Step();
		*datarr=E_grid; ++datarr;
	}
}
// ------------------------------------------------------------------------------------------------
#endif
