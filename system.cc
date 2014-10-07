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
double Cos[N];
double Sin[N];
double Sindot[N];

double T[N][N];
	//double D[N][N];
double T_Cos[N][N];
double T_Sin[N][N];
double T_Sindot[N][N];
double R_Cos[N][N];
double R_Sin[N][N];
double R_Sindot[N][N];


// Diese Funktion laeuft nur einmal, darf deshalb unoekonomisch sein
void TridiagToeplitz() {
	for (int i=0; i<N; i++) {
		// delta_t=h/E_0
		//D[i][i]=((2*k/m)*(cos((i+1)*M_PI/(N+1)) - 1));
		Cos[i]=cos(sqrt(-((2*k/m)*(cos((i+1)*M_PI/(N+1)) - 1)))*h/E_0);// ACHTUNG: assuming matrix D is neg. def.
		   Sin[i]=sin(sqrt(-((2*k/m)*(cos((i+1)*M_PI/(N+1)) - 1)))*h/E_0)/(sqrt(-((2*k/m)*(cos((i+1)*M_PI/(N+1)) - 1)))); // divided by omega
		Sindot[i]=sin(sqrt(-((2*k/m)*(cos((i+1)*M_PI/(N+1)) - 1)))*h/E_0)*(sqrt(-((2*k/m)*(cos((i+1)*M_PI/(N+1)) - 1)))); // multiplied by omega
		
		// Wie kann ich diesen extra loop vermeiden?
		double norm=0;
		for (int j=0; j<N; j++) {
			norm+=sin((j+1)*M_PI*(i+1)/(N+1))*sin((j+1)*M_PI*(i+1)/(N+1));
		}
		norm=sqrt(norm);
		for (int j=0; j<N; j++) {
			T[j][i]=sin((j+1)*M_PI*(i+1)/(N+1))/norm;
			
			// die Spalten der Matrix mit den Faktoren fuer die Zeitentwicklung ergaenzen (also sparsity ausnutzen)
			T_Cos[j][i]=sin((j+1)*M_PI*(i+1)/(N+1))*Cos[i]/norm;
			T_Sin[j][i]=sin((j+1)*M_PI*(i+1)/(N+1))*Sin[i]/norm;
			T_Sindot[j][i]=sin((j+1)*M_PI*(i+1)/(N+1))*Sindot[i]/norm;
		}
	}
	
	// Matrizen berechnen, die Startwerte direkt mit Zeitentwicklung verbinden: z. B. T*Cos*T^-1
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			R_Cos[i][j]=0;
			R_Sin[i][j]=0;
			R_Sindot[i][j]=0;
			for (int k=0; k<N; k++) {
				   R_Cos[i][j]+=   T_Cos[i][k]*T[j][k]; // T transposed...
				   R_Sin[i][j]+=   T_Sin[i][k]*T[j][k]; // T transposed...
				R_Sindot[i][j]+=T_Sindot[i][k]*T[j][k]; // T transposed...
			}
		}
	}
	
}

class System {
   public:
   	// constructor
	System(double, double, double*, double*);
	
	// Teilchen
	double pos;
	double v;
	
	// Gitter
	double x[N];
	double xdot[N];
	
	double w[N];
	double wdot[N];
	
	// pointer
	double* x_ptr;
	double* w_ptr;
	double* xdot_ptr;
	double* wdot_ptr;
	
	// Methoden
	double Collision(double, double, double, double);
	void Oscillate();
	double* Evolve(double*);
};

// initializing double delta_t, double pos, double v, double* x, double* xdot
System::System(double pos_0, double v_0, double* x_0, double* xdot_0) {
	pos=pos_0;
	v=v_0;
	
	for (int i=0; i<N; i++) {
		x[i]=*x_0;
		x_0++;
		
		xdot[i]=*xdot_0;
		xdot_0++;
	}
	
	// pointer initialisieren
	x_ptr=x;
	xdot_ptr=xdot;
	w_ptr=w;
	wdot_ptr=wdot;
}

double System::Collision(double m1, double v1, double m2, double v2) {
	return 2*((m1*v1 + m2*v2)/(m1 + m2)) - v1;
}

void System::Oscillate() {
	for (int i=0; i<N; i++) {
		*w_ptr=0;
		*wdot_ptr=0;
		for (int j=0; j<N; j++) {
			*w_ptr+=(R_Cos[i][j]*(*x_ptr) + R_Sin[i][j]*(*xdot_ptr));
			*wdot_ptr+=(R_Cos[i][j]*(*xdot_ptr) - R_Sindot[i][j]*(*x_ptr));
			x_ptr++;
			xdot_ptr++;
		}
		w_ptr++;
		wdot_ptr++;
		
		x_ptr=x_ptr-N; // set pointers back to start
		xdot_ptr=xdot_ptr-N;
	}
	w_ptr=w_ptr-N; // set pointers back to start
	wdot_ptr=wdot_ptr-N;

	// interchange actual and previous values
	double* help_ptr;
	double* helpdot_ptr;

	help_ptr=x_ptr;
	helpdot_ptr=xdot_ptr;

	x_ptr=w_ptr;
	xdot_ptr=wdot_ptr;
	
	w_ptr=help_ptr;
	wdot_ptr=helpdot_ptr;
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
	v=this->Collision(M, v, m, *(xdot_ptr + index));
	
	*(xdot_ptr + index)=this->Collision(m, *(xdot_ptr + index), M, w); // use w

	// Position updaten
	pos+=v*h/E_0;
	
	// Gitter weiterentwickeln
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
	*arr=b;
	arr++;
	*arr=pos;
	arr++;
	*arr=index;
	arr++;
	*arr=v;
	arr++;
	*arr=*(xdot_ptr + index);
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