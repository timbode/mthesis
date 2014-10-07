#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

// Systemkonstanten
const unsigned int N=10;

// Energielevel
const unsigned int n=1;
// Boxlaenge
const double a=1.0*1e-12;

const double h=6.62606957*1e-34;
const double m_e=9.10938291*1e-31;
const double c=299792458;

// Teilchenmasse
const double M=m_e;

// Energie
const double E_0=n*n*h*h/(8*M*a*a);

// Gitterteilchenmasse
const double m=m_e/N;
const double L=a/(N+1);
const double k=N*N*m*c*c/(a*a);

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
	System(double, double, double*, double*);
	
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
System::System(double pos_0, double v_0, double* x_0, double* xdot_0) {
	delta_t=1.0; // irgendwie initialisieren: als erstes findet ein Stoss statt
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
	
	// neue Geschwindigkeiten berechnen
	double w=v; // temporarily copy v
	v=this->Collision(M, v, m, xdot[index]);
	
	xdot[index]=this->Collision(m, xdot[index], M, w); // use w
	
	// delta_t updaten
	delta_t=1/((M*c*c/h) + (0.5*M*v*v/h));
	
	// Position updaten
	pos+=v*delta_t;
	
	// update delta_t und Position zusammenfassen
	//pos+=2*h/(M*v);
	
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
	double B=(N+1)*L; // B=a...
	int b=floor(pos/B);
	if ((b % 2) == 0) { // gerade
		pos=pos - b*B;
		// v bleibt
	}
	else { // ungerade
		pos=B - (pos - b*B);
		v=-v;
	}
	
	// Position, resultierenden Index, resultierende Geschwindigkeiten des Teilchens und der Gittermasse speichern (v ist dann die Ausgangsgeschw. fuer den folgenden Stoss)
	file1 << b << "   " << pos << "   " << index << "   " << v << "   " << xdot[index] << endl;
	
	// Dringend beachten: Das Teilchen fliegt mit der hier gespeicherten Geschw. v vom vorigen Ort zum hier gespeicherten Ort.
	// Der Index gehoert zum vorigen Ort, weil er vor dem Update der Position gesetzt wird.
	// b gehoert hingegen zur aktuellen, weil upgedateten Position.
	
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

// Anzahl Zeitschritte
int steps=10;

// Anfangsgeschwindigkeit Teilchen
double v_0=n*h/(2*M*a);

// Raeumliche Aufloesung Anfangswerte
const unsigned int resol=2;
double pos_0s[resol]={50.75646,60.75646};//,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};

particle_data << "# N: " << N << endl;
particle_data << "# M: " << M << endl;
particle_data << "# m: " << m << endl;
particle_data << "# v_0: " << v_0 << endl;
particle_data << "# L: " << L << endl; 
particle_data << "# steps: " << steps << endl;  
particle_data << "# resol: " << resol << endl; 
particle_data << "# ----------------------------" << endl;
particle_data << "# ----------------------------" << endl;

for (int ii=0; ii<resol; ii++) {
	// Anfangsposition Teilchen
	double pos_0=pos_0s[ii]*L;

	System sys(pos_0, v_0, x_0, xdot_0);

	for (int i=0; i<steps; i++) {
		sys.Evolve(particle_data, lattice_positions, lattice_velocities);
	}
}

particle_data.close();
lattice_positions.close();
lattice_velocities.close();

return 0;
}
