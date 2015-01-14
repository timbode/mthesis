// Verlet method for a multidimensional grid
#ifndef VERLET_H
#define VERLET_H

#include <vector>
#include <math.h>
#include <omp.h>
#include <algorithm> // max, min

// random numbers
#include <ctime>
#include </home/students/lappet/Documents/Masterarbeit/boost_1_56_0/boost/random.hpp>
#include </home/students/lappet/Documents/Masterarbeit/boost_1_56_0/boost/generator_iterator.hpp>

// namespace
#include "constants.hh"

using namespace std;
using namespace Constants;

//typedef for random number generation; see http://www.boost.org/doc/libs/1_37_0/libs/random/random_demo.cpp
typedef boost::mt19937 base_generator_type;

// ------------------------------------------------------------------------------------------------

class Verlet {
	public:
		Verlet(int, int, double, double); // construction from file
		~Verlet(); // destructor

		int p; // grid number==particle number
		int rep; // repetition number

		double T; // time interval
		double dt; // recursion time

		double* r0; double* r1; double* r2; // positions: three arrays required because Verlet method depends on the last two sets of positions

		double* rdot; // velocities

		// methods
		int Index(int, int, int, int);
		double NearestNeighbours(int, int, int, int);
		double PotentialEnergy(int, int, int, int, int);
		double Step();
		void Evolve();

}; // do not forget the funny semicolon down here...

// init
Verlet::Verlet(int P, int Rep, double T_0, double dt_0) {
	p=P; rep=Rep;
	T=T_0; dt=dt_0;
	r0=new double[Max]; r1=new double[Max]; r2=new double[Max];
	rdot=new double[Max];

	if (rep==0) {

		//----------------------------------------------------------------------------------------------
		// create a random generator for each system
		base_generator_type generator(42u);
		// set seed to system time (probably change this to something else)
		generator.seed(static_cast<unsigned int>(time(0)));

		// Define a uniform random number distribution which produces "double"
		// values between 0 and 1 (0 inclusive, 1 exclusive).
		boost::uniform_real<> uni_dist(-1,1);
		boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
		//-----------------------------------------------------------------------------------------------

		// default initialization: positions must be set
		for (int x=0; x<N_[0]; ++x) {
			for (int y=0; y<N_[1]; ++y) {
				for (int z=min(0, N_[2]); z<max(0, N_[2]); ++z) {
					  r0[this->Index(x, y, z, 0)]=x;   r0[this->Index(x, y, z, 1)]=y;   r0[this->Index(x, y, z, 2)]=z;
					  r1[this->Index(x, y, z, 0)]=x;   r1[this->Index(x, y, z, 1)]=y;   r1[this->Index(x, y, z, 2)]=z;
					  r2[this->Index(x, y, z, 0)]=x;   r2[this->Index(x, y, z, 1)]=y;   r2[this->Index(x, y, z, 2)]=z;
					rdot[this->Index(x, y, z, 0)]=0;   rdot[this->Index(x, y, z, 1)]=0;   rdot[this->Index(x, y, z, 2)]=0;
				}
			}
		}

		// make excitation
		double scale_x=1e-3; // must not be too large - why, exactly?
		double scale_y=1e-3;
		for (int x=1; x<(N_[0]-1); ++x) {
			for (int y=1; y<(N_[1]-1); ++y) {
				for (int z=min(1, N_[2]-1); z<max(1, N_[2]-1); ++z) {
					//cout << x << y << z << '\n';
					r0[this->Index(x, y, z, 0)]+=scale_x*uni();   r0[this->Index(x, y, z, 1)]+=scale_y*uni();   r0[this->Index(x, y, z, 2)]+=0;
					r1[this->Index(x, y, z, 0)]+=scale_x*uni();   r1[this->Index(x, y, z, 1)]+=scale_y*uni();   r1[this->Index(x, y, z, 2)]+=0;
					r2[this->Index(x, y, z, 0)]+=scale_x*uni();   r2[this->Index(x, y, z, 1)]+=scale_y*uni();   r2[this->Index(x, y, z, 2)]+=0;
				}
			}
		}

	}

	else {
		// tricky bug here: while it is true that r2 is overwritten during the first time step,
		// that holds only for the inner components -
		// the boundary components will be randomly initialized if not set below...
		ifstream state_data;
		ostringstream FileNameStream;
		FileNameStream << _DATA_ << "/init/grid_" << p << "_init_chunk_" << rep << ".dat";
		string FileName=FileNameStream.str();
		state_data.open(FileName.c_str());
		int i=0;
		double rr0; double rr1; double rr2;
		while (state_data >> rr0 >> rr1 >> rr2) {
			r0[i]=rr0; r1[i]=rr1; r2[i]=rr2;
			//cout << "construction: " << i << "  " << r0[i] << "  " << r1[i] << '\n';
			++i;
		}

		state_data.close();
	}

}

Verlet::~Verlet() {
	// save current state of the system --- does this state appear twice (before and after)?
	ofstream state_data;
	ostringstream FileNameStream;
	FileNameStream << _DATA_ << "/init/grid_" << p << "_init_chunk_" << rep + 1 << ".dat";
	string FileName=FileNameStream.str();
	state_data.open(FileName.c_str());
	state_data.precision(15); // precision in writing must be high - otherwise there appear discontinuities in the energy when repeating (division by small dt?)
	// write grid positions to file
	for (int i=0; i<Max; ++i) {
		state_data << r0[i] << '\t';
		state_data << r1[i] << '\t';
		state_data << r2[i] << '\n';
		//cout << "destruction: " << i << "  " << r0[i] << "  " << r1[i] << '\n';
	}

	state_data.close();
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
			for (int x=1; x<(N_[0]-1); ++x) { // Postillon: +++ ++x angeblich schneller als x++ +++
				for (int y=1; y<(N_[1]-1); ++y) {
					for (int z=min(1, N_[2]-1); z<max(1, N_[2]-1); ++z) {
						int index=this->Index(x, y, z, alpha);
						r2[index]=2*r1[index] - r0[index] - ((k/m)*dt*dt)*(2*dim*r1[index] - NearestNeighbours(x, y, z, alpha)); // ((k/m)*dt*dt) must be << 1
						rdot[index]=(r2[index] - r1[index])/dt; // (r2[index] - r0[index])/(2*dt)

						// energy
						E+=0.5*m*rdot[index]*rdot[index] + PotentialEnergy(index, x, y, z, alpha);
					}
				}
			}
		}
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
