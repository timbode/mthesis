#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

// Systemkonstanten
const unsigned int N=20;

// Teilchenmasse
const double M=100.0;
// Gitterteilchenmasse
const double m=1.0;

const double L=1.0;
const double k=10000000.0;//const double k=N*N*m*10e16/(((N+1)*L)*((N+1)*L));

// Matrizen
double T[N][N];
double D[N][N];

void TridiagToeplitz() {
	for (int i=0; i<N; i++) {
		D[i][i]=(2*k/m)*(cos((i+1)*M_PI/(N+1)) - 1);
		
		// Wie kann ich diesen extra loop vermeiden?
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
	System(double, double, double, double*, double*);
	
	// Teilchen
	double delta_t;
	
	double pos;
	double v;
	
	// Gitter
	double x[N];
	double xdot[N];
	
	double y[N];
	double ydot[N];
	
	// Methoden
	double Collision(double, double, double, double);
	void Oscillate();
	void Evolve(ofstream&, ofstream&, ofstream&);
};

// initializing double delta_t, double pos, double v, double* x, double* xdot
System::System(double time_step, double pos_0, double v_0, double* x_0, double* xdot_0) {
	delta_t=time_step;
	pos=pos_0;
	v=v_0;
	for (int i=0; i<N; i++) {
		x[i]=*x_0;
		x_0++;
		
		xdot[i]=*xdot_0;
		xdot_0++;
	}
}

double System::Collision(double m1, double v1, double m2, double v2) {
	return 2*((m1*v1 + m2*v2)/(m1 + m2)) - v1;
}

void System::Oscillate() {
	double w1[N]={};
	double w2[N]={};
	for (int i=0; i<N; i++) {
		w1[i]=y[i]*cos(sqrt(-D[i][i])*delta_t) + ydot[i]/(sqrt(-D[i][i]))*sin(sqrt(-D[i][i])*delta_t); // ACHTUNG: assuming matrix D is neg. def.
		w2[i]=ydot[i]*cos(sqrt(-D[i][i])*delta_t) - y[i]*(sqrt(-D[i][i]))*sin(sqrt(-D[i][i])*delta_t);
	}
		
	for (int i=0; i<N; i++) {
		x[i]=0;
		xdot[i]=0;
		for (int j=0; j<N; j++) {
			x[i]+=T[i][j]*w1[j];
			xdot[i]+=T[i][j]*w2[j];
		}
	}
}
// argument files are particle_data, lattice_positions, lattice_velocities
void System::Evolve(ofstream& file1, ofstream& file2, ofstream& file3) {

	// eine Gittermasse herausgreifen
	int index=0;
	index=round(pos/L) - 1;
	
	// Spezialfall Enden
	if (index==N) { 
		index-=1; // rechts
	}
	if (index==-1) {
		index+=1; // links
	}
	
	// Position, resultierenden Index, Ausgangsgeschwindigkeiten des Teilchens und der Gittermasse speichern
	file1 << pos << "   " << index << "   " << v << "   " << xdot[index] << endl;	
	
	// neue Geschwindigkeiten berechnen
	double w=v; // temporarily copy v
	v=this->Collision(M, v, m, xdot[index]);
	
	xdot[index]=this->Collision(m, xdot[index], M, w); // use w
	
	// Position updaten
	pos+=v*delta_t;
	
	// Gitter weiterentwickeln
	for (int i=0; i<N; i++) {
		y[i]=0;
		ydot[i]=0;
		for (int j=0; j<N; j++) {
			y[i]+=T[i][j]*x[j];
			ydot[i]+=T[i][j]*xdot[j];
			
			// Gittermassenpositionen mit abgreifen (vor Evolution)
			if (i == 0) {
				file2 << x[j] << "   ";
			}
			
			// Gittermassengeschwindigkeiten mit abgreifen (nach Stoss, aber vor Evolution)
			if (i == 0) {
				file3 << xdot[j] << "   ";
			}
		}
	}
	
	// Zeilenumbruch in file2 ergaenzen
	file2 << endl;
	
	// Zeilenumbruch in file3 ergaenzen
	file3 << endl;
	
	this->Oscillate(); // sets x and xdot to new values
	
	// Reflektion
	double B=(N+1)*L; // B=D...
	int b=floor(pos/B);
	if ((b % 2) == 0) { // gerade
		pos=pos - b*B;
		// v bleibt
	}
	else { // ungerade
		pos=B - (pos - b*B);
		v=-v;
	}
	
	//cout << pos << " " << v << endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------

int main() {

// Matrizen erzeugen
TridiagToeplitz();

/*
for (int i=0; i<N; i++) {
	for (int j=0; j<N; j++) {
		cout << T[i][j] << " ";
	}
	cout << endl;
}
*/

// Textdateien anlegen und oeffnen
ofstream particle_data, lattice_positions, lattice_velocities;
particle_data.open("particle_data.txt");
lattice_positions.open("lattice_positions.txt");
lattice_velocities.open("lattice_velocities.txt");

// Anfangswerte Gitter
double x_0[N]={};
double xdot_0[N]={};

// "Frequenz"
double time_step=1.0;

// Anfangswerte Teilchen
double pos_0=10.0;
double v_0=420.0;

particle_data << "# N: " << N << endl;
particle_data << "# M: " << M << endl;
particle_data << "# m: " << m << endl;
particle_data << "# v_0: " << v_0 << endl; 
particle_data << "# ----------------------------" << endl;
particle_data << "# ----------------------------" << endl;

System sys(time_step, pos_0, v_0, x_0, xdot_0);

int steps=1000;
for (int i=0; i<steps; i++) {
	sys.Evolve(particle_data, lattice_positions, lattice_velocities);
}

particle_data.close();
lattice_positions.close();
lattice_velocities.close();

return 0;
}
