// namespace for constant vars
#ifndef CONSTANTS_H
#define CONSTANTS_H
namespace Constants {

// ------------------------------------------------------------------------------------------------
// path to data folder - [19:-2] gives folder in exec.sh
const char _DATA_[]="data";
// constants from verlet.hh
const unsigned int dim=2;
// careful: N_[0], N_[1] must be >=3; N_[2] can be either =1 or >=3 (there have to be inner points for the nearest-neighbour concept to work)
const int N_ [3]={101, 101, 1};
const unsigned int Max=3*N_[0]*N_[1]*N_[2]; // length of arrays
const double k=1e2; // spring constant
const double m=1.0; // grid-point mass
const double d=0.1; // grid-point diameter

const double a_X=1.0; const double a_Y=1.0; const double a_Z=0.0;
const char system_type[]="droplet"; // set this to "particle" or "droplet"
// a_X/(N_[0]-1), a_Y/(N_[1]-1) and a_Z/(N_[2]-1) should be the same to ensure homogeneity and isotropy
//const double L=1.0; // for physical lattice
const double L=a_X/(N_[0] - 1); // for approx of membrane (droplet case)

// possible random excitation on grid
const bool excitation=0;
// ------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------
// "double-slit stuff"
const bool wall=0;
const bool stop_if_crashed=1;
const int half_width=10;
const int  left_face_pos=(N_[0]-1)/4 - half_width; // on the x-axis
const int right_face_pos=(N_[0]-1)/4 + half_width;

// slit 1 BELOW the middle line (y==const) and slit 2 ABOVE
const int slit_shift=10;
const int slit_width=20;
const int slit_1_lower=(N_[1]-1)/2 - slit_shift - slit_width; // on the y-axis
const int slit_1_upper=(N_[1]-1)/2 - slit_shift; // note asymmetric definition

const int slit_2_lower=(N_[1]-1)/2 + slit_shift; // on the y-axis
const int slit_2_upper=(N_[1]-1)/2 + slit_shift + slit_width; // note asymmetric definition
// ------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------
// constants from particle.hh and droplet.hh
const double M=250;//N_[0]*N_[0]-4*N_[0]+4; // particle mass
const double D=0.4; // particle diameter
const double f=100; // "bouncing" frequency of particle, so 1/10 for every 10 seconds
// ------------------------------------------------------------------------------------------------

}
#endif
