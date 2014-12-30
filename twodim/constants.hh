// namespace for constant vars
#ifndef CONSTANTS_H
#define CONSTANTS_H
namespace Constants {

// ------------------------------------------------------------------------------------------------
// constants from verlet.hh
const unsigned int dim=2;
// careful: N_[0], N_[1] must be >=3; N_[2] can be either =1 or >=3 (there have to be inner points for the nearest-neighbour concept to work)
const int N_ [3]={10, 10, 1};
const unsigned int Max=3*N_[0]*N_[1]*N_[2]; // length of arrays
const double k=1e3; // spring constant
const double m=1.0; // grid-point mass
const double d=0.1; //  grid-point diameter

const double a_X=1.0;
const double a_Y=1.0;
const double a_Z=0.0;
const char system_type[]="particle"; // set this to "particle" or "droplet"
// a_X/N_[0], a_Y/N_[1] and a_Z/N_[2] should be the same to ensure homogeneity and isotropy
const double L=1.0; // for physical lattice
//const double L=a_X/N_[0]; // for approx of membrane (droplet case)
// ------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------
// constants from particle.hh and droplet.hh
const double M=64.0; // particle mass
const double D=0.4; // particle diameter
const double f=1.0; // "bouncing" frequency of particle, so 1/10 for every 10 seconds
// ------------------------------------------------------------------------------------------------

}
#endif
