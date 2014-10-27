#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <vector>
#include <math.h>

using namespace std;

// Systemkonstanten
const unsigned int N=500; // Anzahl Gittermassen
const unsigned int steps=1000000; // Anzahl Zeitschritte
const unsigned int n=1; // Energielevel
const double a=1.0; // Boxlaenge
const double L=a/(N+1); // Abstand Gittermassen
const double h=1.0; // Wirkungsquantum //6.62606957*1e-34;
const double M=1.0; // Teilchenmasse
const double m=10000000000*M; // Gitterteilchenmasse
const double E_0=n*n*h*h/(8*M*a*a); // Energie
const double k=(N+1)*(N+1)*m*E_0/(2*M*a*a);// Federkonstante // (-1)*(1/((cos(n*M_PI/(N+1)) - 1)))*((m*h*h*pow(M_PI,2)*pow(n,4))/(32*M*M*pow(a,4)))

double v_0=0.0; // Anfangsgeschwindigkeit Teilchen
int start_index=400; // Anfangsposition Teilchen
double excitation=0.00005; // Anfangsanregung Gitter

//------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------

// Vektoren und Matrizen
vector<double> EigVals(N);
vector< vector<double> > T(N,vector<double>(N));

// this function runs only once
void TridiagToeplitz() {
	for (int i=0; i<N; i++) {
		EigVals[i]=(2*k/m)*(cos((i+1)*M_PI/(N+1)) - 1);
		
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
	int index;
	double v;
	
	// Zeit
	double delta_t;
	// aktuelle Geschw.
	double xdot_index;
	
	// Gitter
	double* y;
	double* ydot;
	
	double* w;
	double* wdot;
	
	// Methoden
	double Collision(double, double, double, double);
	double Oscillate(int, double);
	double* Evolve(double*);
};

// initializing
System::System(double index_0, double v_0, double xdot_index_0, double* y_0, double* ydot_0) {
	index=index_0;
	v=v_0;
	
	delta_t=0.0;
	xdot_index=xdot_index_0;
	
	
	y=new double[N];
	ydot=new double[N];
	
	w=new double[N];
	wdot=new double[N];

	y=y_0;
	ydot=ydot_0;	
}

double System::Collision(double m1, double v1, double m2, double v2) {
	return 2*((m1*v1 + m2*v2)/(m1 + m2)) - v1;
	// return 2*v2 - v1;
}

double System::Oscillate(int ind, double del_t) {
	double next=0.0;
	//#pragma omp parallel for reduction(+:next)
	for (int i=0; i<N; i++) {
		w[i]=y[i]*cos(sqrt(-EigVals[i])*del_t) + (ydot[i]/(sqrt(-EigVals[i])))*sin(sqrt(-EigVals[i])*del_t); // Watch out: assuming EigVals are negative
		wdot[i]=ydot[i]*cos(sqrt(-EigVals[i])*del_t) - y[i]*(sqrt(-EigVals[i]))*sin(sqrt(-EigVals[i])*del_t);
		
		// compute next lattice velocity
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
	v=this->Collision(M, v, m, xdot_index);
	
	double ww_pre=xdot_index; // save xdot_index before collision
	xdot_index=this->Collision(m, xdot_index, M, w); // use w
	//double ww_post=xdot_index; // save xdot_index after collision but before oscillation
	
	// ydot updaten
	#pragma omp parallel for
	for (int i=0; i<N; i++) {
		ydot[i]+=T[index][i]*(xdot_index - ww_pre); //Transponierung aber nicht notwendig: T ist sym.
	}

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
	
	// Gitter weiterentwickeln
	xdot_index=this->Oscillate(index, delta_t); // sets x and xdot to new values
	
	// Hier wird der Index gespeichert, an dem der folgende Stoss mit den hiesigen Geschwindigkeiten stattfindet
	*arr=index;
	arr++;
	*arr=v;
	arr++;
	//*arr=ww_pre; // use ww_pre
	//arr++;
	//*arr=ww_post; // use ww_post
	//arr++;
	*arr=xdot_index;
	arr++;
	
	return arr;
	
	// Beachten: Das Teilchen fliegt mit der hier gespeicherten Geschw. v vom vorigen Ort zum hier gespeicherten Ort.
	
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------

int main() {

// Matrizen erzeugen
TridiagToeplitz();

// Raeumliche Aufloesung Anfangswerte
const unsigned int resol=1;
int index_0s[resol]={start_index}; // integers...


// number of elements in each file line
int neefl=3;
vector<double> particle_data_array(neefl*resol*steps, 0.0); // neefl=number of elements in each file line

// Anfangswerte Gitter
double x_0[N]={};
double xdot_0[N]={};
double y_0[N]={};
double ydot_0[N]={};

ydot_0[n-1]=excitation;//0.01;//0.1118;//0.005;

	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			x_0[i]+=T[i][j]*y_0[j];
			xdot_0[i]+=T[i][j]*ydot_0[j];
			
		}
		//cout << xdot_0[i] << "\n";
	}

double energie=0;
	for (int i=0; i<N; i++) {
			energie+=0.5*m*xdot_0[i]*xdot_0[i];
	}
	
cout << "Energie: " << energie <<", E_0: " << E_0 << '\n';
cout << "\n";

// Anfangsgeschwindigkeit erste Gittermasse
double xdot_index_0=0.0;

double* ptr;
ptr=&particle_data_array[0];

for (int i=0; i<resol; i++) {
	// Anfangsposition Teilchen
	int index_0=index_0s[i]; // round(pos_0/L) - 1
	
	for (int j=0; j<N; j++) {
		xdot_index_0+=T[index_0][j]*ydot_0[j];
	}

	System sys(index_0, v_0, xdot_index_0, y_0, ydot_0);

	for (int k=0; k<steps; k++) {
		ptr=sys.Evolve(ptr);
	}
}

// Textdatei anlegen und oeffnen
ofstream particle_data;
particle_data.open("bouncing_particle_data.txt");

particle_data << "# N: " << N << '\n';
particle_data << "# M: " << M << '\n';
particle_data << "# m: " << m << '\n';
particle_data << "# v_0: " << v_0 << '\n';
particle_data << "# L: " << L << '\n'; 
particle_data << "# steps: " << steps << '\n';  
particle_data << "# resol: " << resol << '\n'; 
particle_data << "# ----------------------------" << '\n';
particle_data << "# ----------------------------" << '\n';

// array einlesen
for (int i=0; i<resol*steps; i++) {
	for (int j=i*neefl; j<(i+1)*neefl; j++) { // neefl=number of elements in each file line
		particle_data << particle_data_array[j] << '\t'; // possibly use << setw(15) or so
	}
	particle_data << '\n';
}

particle_data.close();

return 0;
}
