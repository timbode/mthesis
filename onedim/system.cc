#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <omp.h>
#include <math.h>
#include <vector>

#include "system_constants.hh"

using namespace std;
using namespace SystemConstants;

// github update "move and bounce"...

// Vektoren und Matrizen
vector<double> EigVals(N);
vector< vector<double> > T(N,vector<double>(N));

// this function runs only once
void TridiagToeplitz() {
	for (int i=0; i<N; i++) {
		EigVals[i]=(2*k/m)*(cos((i+1)*M_PI/(N+1)) - 1);
		if (i==0) cout << 2*M_PI/sqrt(-EigVals[i]) << '\n';
		double norm=0;
		for (int j=0; j<N; j++) {
			norm+=sin((j+1)*M_PI*(i+1)/(N+1))*sin((j+1)*M_PI*(i+1)/(N+1));
		}
		norm=sqrt(norm);
		for (int j=0; j<N; j++) {
			T[j][i]=sin((j+1)*M_PI*(i+1)/(N+1))/norm;
		}
	}
}

class System {
   public:
   	// constructor
	System(double, double, double, double, double*, double*);
	// destructor
	~System();

	// Teilchen
	double pos;
	int index;
	double v;

	// Zeit
	double delta_t;

	// aktuelle Gittergeschw.
	double xdot_index;

	// Gitter
	double* y;
	double* ydot;

	double* w;
	double* wdot;

	// Methoden
	double Collide(double, double, double, double);
	int Move();
	int Bounce();
	double Oscillate(int, double);
	double* Evolve(double*);
};

// initializing
System::System(double delta_t_0, double pos_0, double v_0, double xdot_index_0, double* y_0, double* ydot_0) {
	pos=pos_0;
	index=round(pos_0/L) - 1;
	v=v_0;

	delta_t=delta_t_0; // (2*k/m)*(cos(M_PI/(N+1)) - 1);
	xdot_index=xdot_index_0;

	y=new double[N];
	ydot=new double[N];

	w=new double[N];
	wdot=new double[N];

	y=y_0;
	ydot=ydot_0;
}

// destructor
System::~System() {

}

double System::Collide(double m1, double v1, double m2, double v2) {
	return 2*((m1*v1 + m2*v2)/(m1 + m2)) - v1;
}

int System::Move() {
	// Position updaten
	pos+=v*delta_t;

	// Reflektion
	int B=floor(pos/a);
	if ((B % 2) == 0) { // gerade
		pos=pos - B*a;
		// v bleibt
	}
	else { // ungerade
		pos=a - (pos - B*a);
		v=-v;
	}

	// Index updaten
	index=round(pos/L) - 1;

	// Spezialfall Enden
	if (index== N) index-=1; // rechts
	if (index==-1) index+=1; // links

	return B;
}

int System::Bounce() {
	// Index updaten
	int sign_v;
	if (v > 0) {
		sign_v=1;
	}
	else if (v < 0) {
		sign_v=-1;
	}
	else {
		sign_v=1e-15;
		cout << "Watch out: v has become zero..." << '\n';
	}
	index+=sign_v;

	delta_t=L/(sign_v*v);

	// Spezialfall Enden
	if (index==N) {
		index-=1; // rechts
		delta_t=2*delta_t;
		v=-v;
	}
	if (index==-1) {
		index+=1; // links
		delta_t=2*delta_t;
		v=-v;
	}

	// update position
	pos=(index + 1)*L; // pos+=sign_v*L

	return 0;
}

double System::Oscillate(int ind, double del_t) {
	double next=0.0;
	#pragma omp parallel for ordered reduction(+:next)
	for (int i=0; i<N; i++) {
		w[i]   =y[i]*cos(sqrt(-EigVals[i])*del_t)+ (ydot[i]/(sqrt(-EigVals[i])))*sin(sqrt(-EigVals[i])*del_t); // Watch out: assuming EigVals are negative
		wdot[i]=ydot[i]*cos(sqrt(-EigVals[i])*del_t) - y[i]*(sqrt(-EigVals[i]))*sin(sqrt(-EigVals[i])*del_t);

		// compute next lattice velocity
		#pragma ordered
		next+=T[ind][i]*wdot[i];
	}

	// interchange up-to-date and previous values
	double* help_ptr;
	double* helpdot_ptr;

	help_ptr=y;
	helpdot_ptr=ydot;

	y=w;
	ydot=wdot;

	w=help_ptr;
	wdot=helpdot_ptr;

	return next;
}

double* System::Evolve(double* arr) {

	// neue Geschwindigkeiten berechnen
	double w=v; // temporarily copy v
	v=this->Collide(M, v, m, xdot_index);

	double ww_pre=xdot_index; // save xdot_index before collision
	xdot_index=this->Collide(m, xdot_index, M, w); // use w
	double ww_post=xdot_index; // save xdot_index after collision but before oscillation

	// ydot updaten
	#pragma omp parallel for
	for (int i=0; i<N; i++) {
		ydot[i]+=T[index][i]*(xdot_index - ww_pre); //Transponierung aber nicht notwendig: T ist sym.
	}

	// move particle
	int b=this->Move();

	// Gitter weiterentwickeln
	xdot_index=this->Oscillate(index, delta_t);

	// Position, resultierenden Index, (resultierende) Geschwindigkeiten des Teilchens und der Gittermasse speichern (v ist dann die Ausgangsgeschw. fuer den folgenden Stoss)
	*arr=b; arr++;
	*arr=pos; arr++;
	*arr=index; arr++;
	*arr=v; arr++; // here, v must be saved for the weighting
	*arr=ww_pre; arr++;// use ww_pre
	*arr=ww_post; arr++;// use ww_post

	return arr;

	// Dringend beachten: Das Teilchen fliegt mit der hier gespeicherten Geschw. v vom vorigen Ort zum hier gespeicherten Ort.
	// Der Index gehoert zum vorigen Ort, weil er vor dem Update der Position gesetzt wird.
	// b gehoert hingegen zur aktuellen, weil upgedateten Position.

}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------

int main() {

// Matrizen erzeugen
TridiagToeplitz();

ofstream system_info;
system_info.open("data/system.txt");

system_info << "# N: " << N << '\n';
system_info << "# M: " << M << '\n';
system_info << "# m: " << m << '\n';
system_info << "# k: " << k << '\n';
system_info << "# L: " << L << '\n';

system_info << "# steps: " << steps << '\n';
system_info << "# chunk_size: " << chunk_size << '\n';
system_info << "# resol: " << resol << '\n';

system_info << "# delta_t: " << delta_t_0 << '\n';
system_info << "# v_0: " << v_0 << '\n';

system_info.close();

// specify starting positions
//double pos_0s[resol];
//for (int i=0; i<resol; i++) pos_0s[i]=i+1;

for (int k=0; k<resol; k++) {
	// must be reset because in the initialization the System pointers are set to these arrays
	double y_0[N]={};
	double ydot_0[N]={};
	ydot_0[n-1]=excitation;

	double pos_0=pos_0s[k]*L; // Anfangsposition Teilchen
	int index_0=round(pos_0/L) - 1; // Anfangsindex Teilchen
	double xdot_index_0=0.0; // Anfangsgeschwindigkeit respektive Gittermasse

	for (int j=0; j<N; j++) {
		xdot_index_0+=T[index_0][j]*ydot_0[j];
	}


	System sys(delta_t_0, pos_0, v_0, xdot_index_0, y_0, ydot_0);

	for (int kk=0; kk<steps/chunk_size; kk++) {

		int neefl=6; // number of elements in each file line
		vector<double> particle_data_array(neefl*chunk_size, 0.0);

		double* ptr;
		ptr=&particle_data_array[0];

		// time evolution
		for (int j=0; j<chunk_size; j++) {
			ptr=sys.Evolve(ptr);
		}

		// open file
		ofstream particle_data;
	    	ostringstream FileNameStream;
	    	FileNameStream << "data/particle_" << k << "_data_" << kk << ".txt";
	    	string FileName=FileNameStream.str();
	    	particle_data.open(FileName.c_str());

		// array einlesen
		for (int i=0; i<chunk_size; i++) {
			for (int j=i*neefl; j<(i+1)*neefl; j++) { // neefl=number of elements in each file line
				particle_data << particle_data_array[j] << '\t'; // possibly use << setw(15) or so
			}
			particle_data << '\n';
		}

		particle_data.close();
	}
}

return 0;
}
