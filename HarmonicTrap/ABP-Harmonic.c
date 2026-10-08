/*
  2026-08-20 
  This code simulates the dynamics of active Brownian particles whose dynamics is given by
  
  \dot x = v_0 \cos\theta - k x + \sqrt{2 D_p} \eta_x 
  \dot y = v_0 \sin\theta - k y  + \sqrt{2 D_p} \eta_y
  \dot \theta = \sqrt{2 D_r} \eta_\theta

  Parameters choice: to test the convergence of the theory as \epsilon=\tau_a/\tau_V \to 0, we set

  The mobility mu=1 so that the potential time scale is \tau_V=1/k.

  We keep the active diffusivity constant,
  
  D_a=1=v_a^2 \tau_a / d = v_r^2 \tau_r / [d (d-1)]
  
  as well as
  
  D_p=1
  
  Our small parameter is eps=\tau_a/\tau_V=k \tau_r/(d-1).
  
  We work in d=2 so that eps=k \tau_r and \tau_r=eps/k and D_r=k/eps.
  
  Once k and \tau_r are fixed, we set v_a from
  
  v_a=\sqrt{d (d-1) D_a/\tau_r}

  The natural length scale of the potential is

  \sqrt{D_a/\mu k}

  To be added, check that
  <r^2> = 2 Dt/k + v_0^2/[ k(k+Dr) ]
  
*/

// This sets the options used in this code
#include "./Options.h"

#define EPS 1e-10
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#ifdef _MT
#include "mt19937-64.c"
#endif

#ifdef _PCG
#include "pcg_variants.h"
#include "pcg_julien.h"
#endif

#include <time.h>

// Structure particles which contains all the data needed to characterize the state of a particle
typedef struct particle{
  double x;
  double y; //position of particles
  double theta; // orientation of particles
} particle;

// Structure param contains the parameters of the code that will be
// passed to functions in a compact way. This starts to be long and
// could be split
typedef struct param {
  long N;  // number of particles
  double dt; // time-step
  double v0; // particle speed
  double Dr; //rotational diffusivity
  double Dt; //translational diffusivity
  double Da; // active diffusivity  
  double sqrt2Drdt; //rotational diffusivity
  double sqrt2Dtdt; //translational diffusivity
  
  // Parameters of the confining potential
  double k; // amplitude of the potential
} param;


/*
  We are computing the number density in 2d, which only depends on the
  radius r, but is normalized according to
  \int d^2r rho(r) = \int_0^\infty 2 \pi r dt \rho(r)=1

  If we call n_i the number of measurements in [r, r+dr], the normalization is

  \sum_i n_i = 1 = 2\pi \int_0^r r dr \rho(r)
  if rho(r) is constant in r_i,r_i+dr, we get for the integral over the ith box
  \rho_i 2\pi(r_{i+1}^2/2-r_i^2). We thus want
  \rho_i = n_i/[N \pi (r_{i+1}^2-r_i^2) ]
*/

typedef struct histo{
  FILE* outputhisto;       // File in which histogram of the density is stored
  double* histogram;       // array in which the density histograms is computed
  double histocount;       // Number of recordings so far
  double dr;               // Histograms is made using bins of width dr
  double rmax;             // Maximal value of rmax used for the histogram
  long Nbin;               // Corresponding number of bins in the histogram
  double NextUpdateHisto;  // Next time at which to store the histogram
  double NextStoreHisto;   // Next time at which to store the histogram  
  double StoreHistoInter;  // Interval between two storage of the histogram
  double UpdateHistoInter; // Interval between two storage of the histogram
} histo;

#include "ABP-Harmonic-functions.c"

int main(int argc, char* argv[]){
  
  /* VARIABLE DECLARATIONe */
  long i;              // counters
  double ell_V;        // Typical potential length
  time_t time_clock;   // current physical time
  double _time;        // current time
  double FinalTime;    // final time of the simulation
  double EquilibTime;  // final time of the simulation  
  particle* Particles; // arrays containing the particles
  param Param;         // Structure containing the parameters
  FILE* outputparam;   // File where parameters are stored  
  histo Histo;         // Structure containing the histogram  
#ifdef _MT
  long long seed;      // seed of the random number generator
#endif
#ifdef _PCG
  pcg128_t seed;       // seed of the random number generator
#endif
  
  /* Read parameters of the simulation*/
  // Take input from command line, check that their number is correct and store them
  CheckAndReadInput(argc,argv, &Param, &seed, &FinalTime, &EquilibTime, &Histo, &outputparam);

#ifdef _MT
  init_genrand64(seed);
#endif
#ifdef _PCG
  pcg64_srandom(seed, (pcg128_t) 1);
#endif
  
  /* Initialize variables */
  _time                  = -1.*EquilibTime;
  Particles              = (particle*) malloc( Param.N * sizeof(particle) );
  
  // Initial conditions
  // Typical potential length
  ell_V           =   sqrt( (Param.Da+Param.Dt) / Param.k); 
  
  for(long i=0; i<Param.N; i++){
    Particles[i].x        = ell_V * (-1.0 + 2 * genrand64_real3() );
    Particles[i].y        = ell_V * (-1.0 + 2 * genrand64_real3() );
    Particles[i].theta    = M_PI  * (-1.0 + 2 * genrand64_real3() );        
  }
  
  printf("Initialisation over\n");
  time_clock       = time(NULL);
  
  /* Run the dynamics */
  
  while(_time<FinalTime){
    
    // Move the particles
    Move_Particles_ABP(Particles,Param);
    
    // Increment time
    _time += Param.dt;
    
    //If it is time, update histogram
    if( _time > Histo.NextUpdateHisto - EPS )
      UpdateHisto( Param , &Histo, Particles );
    
    //If it is time, store the histogram
    if( _time > Histo.NextStoreHisto - EPS )
      StoreHisto( Param, &Histo, _time);
  }
  
  printf("Simulation time: %ld seconds\n",time(NULL)-time_clock);
  fprintf(outputparam,"#Simulation time: %ld seconds\n",time(NULL)-time_clock);
  
  free(Histo.histogram);
  free(Particles);
  fclose(outputparam);
  fclose(Histo.outputhisto);

  return 0;
}


