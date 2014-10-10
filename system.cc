#include <iostream>
#include <fstream>
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
const double E_0=n*n*h*h/(8*M*a*a);

// Gitterteilchenmasse
// ca. 1000*M < m < 10000*M bei N=300 und n=1
const double m=10000*M;//1600*M;

const double L=a/(N+1);
const double k=(-1)*(1/((cos(n*M_PI/(N+1)) - 1)))*((m*h*h*pow(M_PI,2)*pow(n,4))/(32*M*M*pow(a,4)));

// Matrizen
double Cos[N];
double Sin[N];
double Sindot[N];

double T[N][N];
double T_Cos[N][N];
double T_Sin[N][N];
double T_Sindot[N][N];

double R_Cos[N][N];
double R_Sin[N][N];
double R_Sindot[N][N];

// Diese Funktion laeuft nur einmal, darf deshalb unoekonomisch sein
void TridiagToeplitz() {
	cout << k << '\n';
	cout << h/E_0 << '\n';
	cout << 2*M_PI << "   " << sqrt(-((2*k/m)*(cos(n*M_PI/(N+1)) - 1)))/E_0<< '\n';
	for (int i=0; i<N; i++) {
		double delta_t=h/E_0;
		double EigVal=((2*k/m)*(cos((i+1)*M_PI/(N+1)) - 1));
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
double pos_0s[resol]={50.0};//,60.75646,70.75646,80.75646,99.75646};

// Anzahl Zeitschritte
const unsigned int steps=100000;
vector<double> particle_data_array(5*resol*steps, 0.0);

// Anfangswerte Gitter
double x_0[N]={};
double xdot_0[N]={};
double y_0[N]={};
double ydot_0[N]={};

ydot_0[n-1]=0.01;//0.005;

	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			x_0[i]+=T[i][j]*y_0[j];
			xdot_0[i]+=T[i][j]*ydot_0[j];
		}
	}

double energie=0;
	for (int i=0; i<N; i++) {
			energie+=0.5*m*xdot_0[i]*xdot_0[i];
	}
	
cout << "Energie: " << energie <<", E_0: " << E_0 << '\n';

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
	for (int j=i*5; j<(i+1)*5; j++) {
		particle_data << particle_data_array[j] << "   ";
	}
	particle_data << '\n';
}

particle_data.close();

return 0;
}
