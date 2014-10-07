#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

// Systemkonstanten
const unsigned int N=10;

// Energielevel
const unsigned int n=1;
// Boxlaenge
const double a=1.0;

const double h=1.0;//6.62606957*1e-34;
const double m_e=9.10938291*1e-31;
const double c=299792458;

// Teilchenmasse
const double M=1.0;

// Energie
const double E_0=1000.0;//n*n*h*h/(8*M*a*a);

// Gitterteilchenmasse
const double m=M;
const double L=a/(N+1);
const double k=1000;

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
	double* Evolve(double*);
};

// initializing double delta_t, double pos, double v, double* x, double* xdot
System::System(double pos_0, double v_0, double* x_0, double* xdot_0) {
	delta_t=h/E_0; // irgendwie initialisieren: als erstes findet ein Stoss statt
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
	for (int i=0; i<N; i++) {
		x[i]=0;
		xdot[i]=0;
		for (int j=0; j<N; j++) {
			x[i]+=T[i][j]*(y[j]*cos(sqrt(-D[j][j])*delta_t) + (ydot[j]/(sqrt(-D[j][j])))*sin(sqrt(-D[j][j])*delta_t)); // ACHTUNG: assuming matrix D is neg. def.
			xdot[i]+=T[i][j]*(ydot[j]*cos(sqrt(-D[j][j])*delta_t) - y[j]*(sqrt(-D[j][j]))*sin(sqrt(-D[j][j])*delta_t));
		}
	}
}

// argument files are particle_data, lattice_positions, lattice_velocities
double* System::Evolve(double* arr) {

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
	//delta_t=1/((M*c*c/h) + (0.5*M*v*v/h));
	
	// Position updaten
	pos+=v*delta_t;
	
	// update delta_t und Position zusammenfassen
	//pos+=2*h/(M*v);
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------
	
	// Gitter weiterentwickeln
	for (int i=0; i<N; i++) {
		y[i]=0;
		ydot[i]=0;
		for (int j=0; j<N; j++) {
			y[i]+=T[j][i]*x[j]; // Transponierung aber nicht notwendig: T ist sym.
			ydot[i]+=T[j][i]*xdot[j];
		}
	}
	
	this->Oscillate(); // sets x and xdot to new values
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------
	
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
	//file1 << b << "   " << pos << "   " << index << "   " << v << "   " << xdot[index] << endl;
	
	*arr=b;
	arr++;
	*arr=pos;
	arr++;
	*arr=index;
	arr++;
	*arr=v;
	arr++;
	*arr=xdot[index];
	arr++;
	
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

/*
for (int i=0; i<N; i++) {
	for (int j=0; j<N; j++) {
		cout << T[i][j] << " ";
	}
	cout << endl;
}
*/

// Raeumliche Aufloesung Anfangswerte
const unsigned int resol=1;
double pos_0s[resol]={9.0};//,60.75646,70.75646,80.75646,90.75646};

// Anzahl Zeitschritte
const unsigned int steps=10;
double particle_data_array[5*resol*steps]={}; // Passt das?

// Anfangswerte Gitter
double x_0[N]={};
double xdot_0[N]={};

// Anfangsgeschwindigkeit Teilchen
double v_0=10000.0;

double* ptr;
ptr=particle_data_array;

for (int ii=0; ii<resol; ii++) {
	// Anfangsposition Teilchen
	double pos_0=pos_0s[ii]*L;

	System sys(pos_0, v_0, x_0, xdot_0);

	for (int i=0; i<steps; i++) {
		ptr=sys.Evolve(ptr);
	}
}

// Textdatei anlegen und oeffnen
ofstream particle_data;
particle_data.open("particle_data.txt");

particle_data << "# N: " << N << endl;
particle_data << "# M: " << M << endl;
particle_data << "# m: " << m << endl;
particle_data << "# v_0: " << v_0 << endl;
particle_data << "# L: " << L << endl; 
particle_data << "# steps: " << steps << endl;  
particle_data << "# resol: " << resol << endl; 
particle_data << "# ----------------------------" << endl;
particle_data << "# ----------------------------" << endl;

// array einlesen
for (int i=0; i<resol*steps; i++) {
	for (int j=i*5; j<(i+1)*5; j++) {
		particle_data << particle_data_array[j] << "   ";
	}
	particle_data << endl;
}

particle_data.close();



return 0;
}
