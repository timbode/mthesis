// namespace for constant vars
#ifndef CONSTANTS_H
#define CONSTANTS_H
namespace Constants {

// ------------------------------------------------------------------------------------------------
// constants from verlet.hh
const unsigned int dim=2;
// careful: N_[0], N_[1] must be >=3; N_[2] can be either =1 or >=3 (there have to be inner points for the nearest-neighbour concept to work properly)
const int N_ [3]={100, 100, 1};
const unsigned int Max=3*N_[0]*N_[1]*N_[2]; // length of arrays
const double k=1e0; // spring constant
const double m=10.0; // grid-point mass
const double d=0.1; //  grid-point diameter
// ------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------
// constants from particle.hh
const double M=100.0; // particle mass
const double D=0.4; // particle diameter
// ------------------------------------------------------------------------------------------------

}
#endif
