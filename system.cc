#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

// Systemkonstanten
const unsigned int N=10;

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

// Vektoren und Matrizen
boost::numeric::ublas::vector<double> Cos (N);
boost::numeric::ublas::vector<double> Sin (N);
boost::numeric::ublas::vector<double> Sindot (N);

boost::numeric::ublas::matrix<double> T (N, N);
boost::numeric::ublas::matrix<double> T_Cos (N, N);
boost::numeric::ublas::matrix<double> T_Sin (N, N);
boost::numeric::ublas::matrix<double> T_Sindot (N, N);

boost::numeric::ublas::matrix<double> R_Cos (N, N);
boost::numeric::ublas::matrix<double> R_Sin (N, N);
boost::numeric::ublas::matrix<double> R_Sindot (N, N);


// Diese Funktion laeuft nur einmal
void TridiagToeplitz() {
	cout << setprecision(10) << k << '\n';
	cout << h/E_0 << '\n';
	cout << 2*M_PI << "   " << sqrt(-((2*k/m)*(cos(n*M_PI/(N+1)) - 1)))/E_0<< '\n';
	for (int i=0; i<N; i++) {
		double delta_t=h/E_0;
		double EigVal=((2*k/m)*(cos((i+1)*M_PI/(N+1)) - 1));
		Cos(i)=cos(sqrt(-EigVal)*delta_t);// ACHTUNG: assuming matrix D is neg. def.
		   Sin(i)=sin(sqrt(-EigVal)*delta_t)/(sqrt(-EigVal)); // divided by omega
		Sindot(i)=sin(sqrt(-EigVal)*delta_t)*(sqrt(-EigVal)); // multiplied by omega
		
		// Wie kann ich diesen extra loop vermeiden?
		double norm=0;
		for (int j=0; j<N; j++) {
			norm+=sin((j+1)*M_PI*(i+1)/(N+1))*sin((j+1)*M_PI*(i+1)/(N+1));
		}
		norm=sqrt(norm);
		for (int j=0; j<N; j++) {
			T(j,i)=sin((j+1)*M_PI*(i+1)/(N+1))/norm;
			
			// die Spalten der Matrix mit den Faktoren fuer die Zeitentwicklung ergaenzen (also sparsity ausnutzen)
			T_Cos(j,i)=T(j,i)*Cos(i);
			T_Sin(j,i)=T(j,i)*Sin(i);
			T_Sindot(j,i)=T(j,i)*Sindot(i);
		}
	}
	
	
	// Matrizen berechnen, die Startwerte direkt mit Zeitentwicklung verbinden: z. B. T*Cos*T^-1
	R_Cos=boost::numeric::ublas::prod(T_Cos,trans(T)); // trans could be omitted (T is sym.)
	R_Sin=boost::numeric::ublas::prod(T_Sin,trans(T));
	R_Sindot=boost::numeric::ublas::prod(T_Sindot,trans(T));
	
}

class System {
   public:
   	// constructor
	System(double, double, boost::numeric::ublas::vector<double>, boost::numeric::ublas::vector<double>);
	
	// Teilchen
	double pos;
	double v;
	
	// Gitter
	boost::numeric::ublas::vector<double> x;
	boost::numeric::ublas::vector<double> xdot;
	
	boost::numeric::ublas::vector<double> w;
	boost::numeric::ublas::vector<double> wdot;
	
	// Methoden
	double Collision(double, double, double, double);
	void Oscillate();
	double* Evolve(double*);
};

// initializing
System::System(double pos_0, double v_0, boost::numeric::ublas::vector<double> x_0, boost::numeric::ublas::vector<double> xdot_0) {
	pos=pos_0;
	v=v_0;

	x=x_0;
	xdot=xdot_0;
}

double System::Collision(double m1, double v1, double m2, double v2) {
	return 2*((m1*v1 + m2*v2)/(m1 + m2)) - v1;
}

void System::Oscillate() {
	w=boost::numeric::ublas::prod(R_Cos,x) + boost::numeric::ublas::prod(R_Sin,xdot);
	wdot=boost::numeric::ublas::prod(R_Cos,xdot) - boost::numeric::ublas::prod(R_Sindot,x);
	x=w;
	xdot=wdot;
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
	v=this->Collision(M, v, m, xdot(index));
	
	xdot(index)=this->Collision(m, xdot(index), M, w); // use w

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
	*arr=xdot(index);
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
const unsigned int steps=1000;
vector<double> particle_data_array(5*resol*steps, 0.0);

// Anfangswerte Gitter
boost::numeric::ublas::vector<double> x_0 (N);
boost::numeric::ublas::vector<double> xdot_0 (N);
boost::numeric::ublas::vector<double> y_0 (N);
boost::numeric::ublas::vector<double> ydot_0 (N);

ydot_0(n-1)=0.1;//0.005;
x_0=boost::numeric::ublas::prod(T,y_0);
xdot_0=boost::numeric::ublas::prod(T,ydot_0);

double energie=0;
	for (int i=0; i<N; i++) {
			energie+=0.5*m*xdot_0(i)*xdot_0(i);
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
