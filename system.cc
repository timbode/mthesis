#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <vector>
#include <math.h>

using namespace std;

// Systemkonstanten
const unsigned int N=500;

// Energielevel
const unsigned int n=1;

// Boxlaenge
const double a=1.0;

// Wirkungsquantum
const double h=1.0;//6.62606957*1e-34;

// Teilchenmasse
const double M=1.0;//0.125;

// Energie
const double E_0=n*n*h*h/(8*M*a*a);//0.125

// Gitterteilchenmasse
const double m=10000*M;//100000000*M;

// Anzahl Zeitschritte
const unsigned int steps=1000000;

// Wiederkehrdauer
const double delta_t=34.0;//506.0;// mind. t=16+2=18

const double L=a/(N+1);
const double k=(N+1)*(N+1)*m*E_0/(2*M*a*a);//(-1)*(1/((cos(n*M_PI/(N+1)) - 1)))*((m*h*h*pow(M_PI,2)*pow(n,4))/(32*M*M*pow(a,4)))

// Matrizen
vector<double>    Cos(N);
vector<double>    Sin(N);
vector<double> Sindot(N);

vector< vector<double> >        T(N,vector<double>(N));
vector< vector<double> >    T_Cos(N,vector<double>(N));
vector< vector<double> >    T_Sin(N,vector<double>(N));
vector< vector<double> > T_Sindot(N,vector<double>(N));

vector< vector<double> >    R_Cos(N,vector<double>(N));
vector< vector<double> >    R_Sin(N,vector<double>(N));
vector< vector<double> > R_Sindot(N,vector<double>(N));

// Diese Funktion laeuft nur einmal
void TridiagToeplitz() {
	cout << "------------------------------------------------\n";
	cout << setprecision(10) << k << '\n';
	cout << h/E_0 << '\n';
	cout << 2*M_PI << "   " << sqrt(-((2*k/m)*(cos(n*M_PI/(N+1)) - 1)))/E_0 << "   " << (sqrt(-((2*k/m)*(cos(n*M_PI/(N+1)) - 1)))*delta_t)/(2*M_PI) << '\n';
	for (int i=0; i<N; i++) {
		//double delta_t=1.0;//h/E_0;
		double EigVal=((2*k/m)*(cos((i+1)*M_PI/(N+1)) - 1));
		//cout << i << "th " << "Eigenvalue: " << sqrt(-EigVal) << '\n';
		Cos[i]=cos(sqrt(-EigVal)*delta_t);// ACHTUNG: assuming matrix D is neg. def.
		   Sin[i]=sin(sqrt(-EigVal)*delta_t)/(sqrt(-EigVal)); // divided by omega
		Sindot[i]=sin(sqrt(-EigVal)*delta_t)*(sqrt(-EigVal)); // multiplied by omega
		
		// Wie kann ich diesen extra loop vermeiden?
		double norm=0;
		for (int j=0; j<N; j++) {
			norm+=sin((j+1)*M_PI*(i+1)/(N+1))*sin((j+1)*M_PI*(i+1)/(N+1));
		}
		norm=sqrt(norm);
		for (int j=0; j<N; j++) {
			T[j][i]=sin((j+1)*M_PI*(i+1)/(N+1))/norm;
			
			// die Spalten der Matrix mit den Faktoren fuer die Zeitentwicklung ergaenzen (also sparsity ausnutzen)
			T_Cos[j][i]=T[j][i]*Cos[i];
			T_Sin[j][i]=T[j][i]*Sin[i];
			T_Sindot[j][i]=T[j][i]*Sindot[i];
		}
	}
	
	// Matrizen berechnen, die Startwerte direkt mit Zeitentwicklung verbinden: z. B. T*Cos*T^-1
	//parallel
	#pragma omp parallel for
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
	double* x;
	double* xdot;
	
	double* w;
	double* wdot;
	
	// Methoden
	double Collision(double, double, double, double);
	void Oscillate();
	double* Evolve(double*);
};

// initializing double delta_t, double pos, double v, double* x, double* xdot
System::System(double pos_0, double v_0, double* x_0, double* xdot_0) {
	pos=pos_0;
	v=v_0;
	
	x=new double[N];
	xdot=new double[N];
	
	w=new double[N];
	wdot=new double[N];

	x=x_0;
	xdot=xdot_0;	
}

double System::Collision(double m1, double v1, double m2, double v2) {
	return 2*((m1*v1 + m2*v2)/(m1 + m2)) - v1;
	// return 2*v2 - v1;
}

void System::Oscillate() {
	#pragma omp parallel for
	for (int i=0; i<N; i++) {
		w[i]=0.0;
		wdot[i]=0.0;
		for (int j=0; j<N; j++) {
			w[i]+=(R_Cos[i][j]*x[j] + R_Sin[i][j]*xdot[j]);
			wdot[i]+=(R_Cos[i][j]*xdot[j] - R_Sindot[i][j]*x[j]);
		}
	}
	// interchange actual and previous values
	double* help_ptr;
	double* helpdot_ptr;

	help_ptr=x;
	helpdot_ptr=xdot;

	x=w;
	xdot=wdot;
	
	w=help_ptr;
	wdot=helpdot_ptr;
}

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
	
	double ww_pre=xdot[index]; // save xdot[index] before collision
	xdot[index]=this->Collision(m, xdot[index], M, w); // use w
	double ww_post=xdot[index]; // save xdot[index] after collision but before oscillation

	// Position updaten
	pos+=v*delta_t;//h/E_0;
	
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
	
	// Position, resultierenden Index, (resultierende) Geschwindigkeiten des Teilchens und der Gittermasse speichern (v ist dann die Ausgangsgeschw. fuer den folgenden Stoss)
	*arr=b;
	arr++;
	*arr=pos;
	arr++;
	*arr=index;
	arr++;
	*arr=v;
	arr++;
	*arr=ww_pre; // use ww_pre
	arr++;
	*arr=ww_post; // use ww_post
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
double pos_0s[resol]={36.0};

vector<double> particle_data_array(6*resol*steps, 0.0); // 6=number of elements in each file line

// Anfangswerte Gitter
double x_0[N]={};
double xdot_0[N]={};
double y_0[N]={};
double ydot_0[N]={};

ydot_0[n-1]=0.005;//15.0;//0.1118;//0.005;

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

// Anfangsgeschwindigkeit Teilchen
double v_0=0.0;//sqrt(2*E_0/M);

double* ptr;
ptr=&particle_data_array[0];

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
	for (int j=i*6; j<(i+1)*6; j++) { // 6=number of elements in each file line
		particle_data << particle_data_array[j] << '\t'; // possibly use << setw(15) or so
	}
	particle_data << '\n';
}

particle_data.close();

return 0;
}
