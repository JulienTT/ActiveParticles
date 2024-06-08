//#define _MT
#define _PCG

#define HARMONIC 
/*
  We use harmonic spheres F = -\grad V where
  
  V(r)= k (r-sigma)^2

  and 

  F_i(r_i-r_j) = 2k (sigma/r_ij-1) (r_i-r_j)
  
  We use a cut-off at r=sigma
  
  There is for now no adaptative time-stepping because of the order dt and \sqrt{dt}.  
*/

//#define LJ 
/*
  We use Lennard-Jones interactions for F = -\grad V where
  
  V=4 epsilon (sigma^12/r^12-sigma^6/r^6)
  
  We use a cut-off at r=sigma*2.7
  
  There is for now no adaptative time-stepping because of the order dt and \sqrt{dt}.  
*/

//#define WCA
/*
  We use WCA potential for F = -\grad V where
  
  V=4 epsilon (sigma^12/r^12-sigma^6/r^6)
  
  We use the LJ potential with a cut-off at r=sigma*2^1/6
  
  There is for now no adaptative time-stepping because of the order dt and \sqrt{dt}.
*/

/*
  Choose the initial condition (IC) between random IC, if RANDOMIC is
  defined, or read from an input if GIVENIC is defined, or SLABIC to
  start in a slab configuration
*/

#define RANDOMIC
//#define GIVENIC
//#define SLABIC

/*
  Choose the boundary conditions between periodic, if PBC is defined,
  and closed along x, if CLOSEDBC is defined
  
  If closed boundary conditions are used, then a confining potential
  is inserted INSIDE the box, to keep the particles far away from the
  periodic boundary condition of the cell. If particles get close
  enough to the wall, the program protests.
  
*/
#define PBC
//#define CLOSEDBC

/*
  Uncomment to compute the pressure exerted on right and left walls
*/
//#define PRESSURE

/*
  Uncomment to compute average of some observables. So far we compute
  the average kinetic energy
*/
//#define OBSERVABLE

/*
  Uncomment to compute average of the stress tensor
*/
#define STRESSTENSOR

/*
  Turn on if the speed is not homogeneous.
  The function Vofr should then be defined
*/
#define VOFR
