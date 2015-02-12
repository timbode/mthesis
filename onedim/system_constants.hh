#ifndef SYSTEM_CONSTANTS_H
#define SYSTEM_CONSTANTS_H
namespace SystemConstants {

// github update "move and bounce"...

// Systemkonstanten
const unsigned int N=300; // Anzahl Gittermassen
const unsigned int steps=5e5; // Anzahl Zeitschritte
const unsigned int chunk_size=steps/10; // size of chunks the evolution is broken-up into
const unsigned int n=1; // Energielevel
const double a=1.0; // Boxlaenge
const double L=a/(N+1); // Abstand Gittermassen
const double h=1.0; // Wirkungsquantum //6.62606957*1e-34;
const double M=1.0; // Teilchenmasse
const double m=1e+4*M; // Gitterteilchenmasse
const double E_0=n*n*h*h/(8*M*a*a); // Energie
const double k=(-1)*(1/((cos(n*M_PI/(N+1)) - 1)))*((m*h*h*pow(M_PI,2)*pow(n,4))/(32*M*M*pow(a,4)));//(N+1)*(N+1)*m*E_0/(2*M*a*a);// Federkonstante

double v_0=0.0; // Anfangsgeschwindigkeit Teilchen
int start_index=36.0; // Anfangsposition Teilchen
const unsigned int resol=1; // Raeumliche Aufloesung Anfangswerte
double pos_0s[resol]={start_index};
double excitation=0.5; // Anfangsanregung Gitter

double delta_t_0=6.0; // mind. 16+2=18 --- stimmt nicht: 6 mit Energie 0.1 funktioniert auch

}
#endif
