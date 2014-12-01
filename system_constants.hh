#ifndef SYSTEM_CONSTANTS_H
#define SYSTEM_CONSTANTS_H
namespace SystemConstants {

const unsigned int N=1000; // Anzahl Gittermassen
const double M=1.0; // Teilchenmasse
const double m=50*M/N; // Gitterteilchenmasse
const double a=1.0; // box length
const double k=1e3;// Federkonstante
const double L=a/(N+1); // Abstand Gittermassen

const unsigned int steps=1e6; // Anzahl Zeitschritte
const unsigned int chunk_size=steps/1; // size of chunks the evolution is broken-up into
const unsigned int resol=N; // Raeumliche Aufloesung Anfangswerte

const double delta_t_0=0.0; // Wiederkehrdauer

const double v_0=0.001; // Anfangsgeschwindigkeit Teilchen
const double excitation=0.0; // Anfangsanregung Gitter


}
#endif
