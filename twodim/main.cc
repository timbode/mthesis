// main file
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>

#include "verlet.hh"
#include "particle.hh"
#include "droplet.hh"

using namespace std;
using namespace Constants;

unsigned long long rdtsc(){
	unsigned int lo,hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return ((unsigned long long)hi << 32) | lo;
}

int main(int argc, char* argv[]) {

int repeat=atoi(argv[1]);
int p=atoi(argv[2]);
int rep=atoi(argv[3]);

int steps=5e5;
double dt=1e-4; // should be 1e-5 or 1e-6
double T=steps*dt;

//----------------------------------------------------------------------------------------------
// create a random generator for each system
base_generator_type GEN(42u);
// set seed to system time (probably change this to something else)
GEN.seed(static_cast<unsigned int>(rdtsc()));

// Define a uniform random number distribution which produces "double"
// values between 0 and 1 (0 inclusive, 1 exclusive).
boost::uniform_real<> UNI_DIST(-1,1);
boost::variate_generator<base_generator_type&, boost::uniform_real<> > UNI(GEN, UNI_DIST);
//-----------------------------------------------------------------------------------------------

//double R_0 [3]={UNI(), UNI(), 0.0}; // watch out: the vectors here MUST NOT be "perfect" (because of the cross product)
//double R_0 [3]={0.01, 0.7 + 1e-1*UNI(), 0.0};
double R_0 [3]={0.01, 0.5 + 2e-1*UNI(), 0.0};
//double R_0 [3]={0.5, 0.52, 0.0};
if (rep==0) {
	cout << "========================================================" << '\n';
	cout << "Start position: ";
	for (int i=0; i<3; ++i) cout << R_0[i] << '\t';
	cout << '\n' << '\n';
	if (wall) {
		cout << "Wall is " << L*2*half_width << " thick" << '\n';
		cout << "Slit centers at " << L*slit_1_lower + (L*slit_width)*0.5 << " (slit 1) and " << L*slit_2_upper - (L*slit_width)*0.5 << " (slit 2)" << '\n';
		cout << "Slit width is " << L*slit_width << '\n';
	}
	cout << "========================================================" << '\n';
}
double V_0 [3]={1.0, 0.0, 0.0};
//double V_0 [3]={0, 0, 0};

unsigned int stats=3;

ofstream system_data;
ostringstream FileNameStream;
FileNameStream << _DATA_ << "/system.dat";
string FileName=FileNameStream.str();
system_data.open(FileName.c_str());
{
	system_data << "# folder: " << _DATA_ << '\n';
	system_data << '\n';
	system_data << "# dim: " << dim << '\n';
	system_data << "# N_X: " << N_[0] << '\n';
	system_data << "# N_Y: " << N_[1] << '\n';
	system_data << "# N_Z: " << N_[2] << '\n';
	system_data << '\n';
	system_data << "# k: " << k << '\n';
	system_data << "# m: " << m << '\n';
	system_data << "# d: " << d << '\n';
	system_data << "# L: " << L << '\n';
	system_data << '\n';
	system_data << "# M: " << M << '\n';
	system_data << "# D: " << D << '\n';
	system_data << "# f: " << f << '\n';
	system_data << '\n';
	system_data << "# steps: " << steps << '\n';
	system_data << "# repeat: " << repeat << '\n';
	system_data << "# dt: " << dt << '\n';
	system_data << "# stats: " << stats << '\n';
	system_data << '\n';
	system_data << "# system_type: " << system_type << '\n';
}
system_data.close();

// keep track of repetition status
cout << " ... " << rep;

// create grid instance
Verlet grid(p, rep, T, dt); // replace T by steps

// make burn-in
if (wall && excitation) {
	grid.Burn_in();
}

vector<double> data_array((dim+2)*steps, 0.0);

if (system_type[0]=='p') {
	// create particle instance
	Particle particle(p, rep, steps, dt, R_0, V_0);
	// give grid to particle and go
	particle.Evolve(&grid, &data_array[0]);
}
else if (system_type[0]=='d') {
	// create droplet instance
	Droplet droplet(p, rep, steps, dt, R_0, V_0);
	// give grid to droplet and go
	droplet.Evolve(&grid, &data_array[0]);
}

// open file
ofstream data;
ostringstream FileNameStream2;
FileNameStream2 << _DATA_ << "/chunks/" << system_type << "_" << p << "_chunk_" << rep << ".dat";
string FileName2=FileNameStream2.str();
data.open(FileName2.c_str());
//data.precision(15);

// write array to file
for (int t=0; t<steps; ++t) {
	for (unsigned int u=0; u<(dim+2); ++u) {
		data << data_array[(dim+2)*t+u] << '\t';
	}
	data << '\n';
}

data.close();

return 0;
}
