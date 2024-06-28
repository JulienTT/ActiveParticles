/*
  2021-02-24 
  This code simulates the dynamics of active Brownian particles whose dynamics is given by
  
  \dot x = v_0 \cos\theta - \partial_x F + \sqrt{2 D_t} \eta_x 
  \dot y = v_0 \sin\theta - \partial_y F + \sqrt{2 D_t} \eta_y
  \dot \theta = \sqrt{2 D_r} \eta_\theta
*/

// This sets the options used in this code
#include "./Options.h"

#define EPS 1e-10
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
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
  long bi; // position of box along x axis
  long bj; // position of box along y axis
} particle;

// Structure used to compute forces in the presence of periodic boundary conditions
typedef struct box{
  long i; // position of box along x axis
  long j; // position of box along y axis
  double epsilonx; //distance to add along x to take into account PBC
  double epsilony; //distance to add along y to take into account PBC
} box;

// Structure param contains the parameters of the code that will be
// passed to functions in a compact way. This starts to be long and
// could be cut
typedef struct param {
  long N;  // number of particles
  double Lx; // system size
  double Ly;
  double dt; // time-step
  double v0; // particle speed
  double Dr; //rotational diffusivity
  double Dt; //translational diffusivity
  double sqrt2Drdt; //rotational diffusivity
  double sqrt2Dtdt; //translational diffusivity

  // Parameters of the interaction potential
  double epsilon; // amplitude of the potential
  double sigma; //translational diffusivity
  double sigma2; //sigma^2
#if defined(WCA) || defined(LJ)
  double sigma6; //sigma^6
  double rmax2;  // interaction cut-off (squared)
  double amp; //used to compute the force
#endif
#if defined(HARMONIC)
  double sigma; //sigma
  double k;  // stiffness
  double ksigma; //used to compute the force
#endif
  
  double mu; //particle mobility
  
  // Parameters of the spatial hashing
  double rbox; //Box size used for hashing. Equal rmax, a priori.
  double rbox2;
  long NxBox; // number of boxes along x
  long NyBox; // number of boxes along y

#ifdef CLOSEDBC
  /*
    The confining potential is Omega/nu (x_wall_left-x)^nu and (x-x_wall_right)^nu
  */
  double x_wall_left;
  double x_wall_right;
  double Omega;
  double nu;
#endif
} param;

#include "ABP-PFAPS-QS-functions.c"

int main(int argc, char* argv[]){
  
  /* VARIABLE DECLARATIONe */
  long i,j,k,i1,j1;    // counters
  time_t time_clock;   // current physical time
  double _time;        // current time
  double rho0;         // average density in the system
  char command_base[1000]=""; // string that contains the desired format of the command line
  
  double FinalTime;    // final time of the simulation
  particle* Particles; // arrays containing the particles
  param Param;  // Structure containing the parameters
#ifdef _MT
  long long seed;      // see of the random number generator
#endif
#ifdef _PCG
  pcg128_t seed;       // see of the random number generator
#endif
  
  double* forces ;     // forces[2*i] is the force along x on part. i
	               // forces[2*i+1] is the force along y on
	               // particle i.
  
  long** Box;          // array containing the first particles in the
                       // box Box[i][j]=k means the first particle in
                       // box i,j is Particle[k]
  
  long* Neighbours;    // array containing the neighbours of each
		       // particle in a box. Neighbours[2*i+1]=k means
		       // the next particle after i is k. k=-1 means i
		       // is the last particle. Neighbours[2*i]=k
		       // means the particle before i is k. k=-1 means
		       // the particle is the first in the box.

  box *** NeighbouringBoxes; //Neighbouringboxes[i][j] contains the
			     //list of neighbouring boxes of box (i,j)

  /* SPECIFIC TO HISTOGRAMS
     To construct a histogram of the data, we average over a period
     StoreHistoInter in an array histogram. We bin together several
     neighbouring boxes, over a linear size width, so that the
     resolution is 1/width^2.
     
     At the end, we plot rho*sigma^2 so that the resolution is  (sigma/width)^2
  */
  FILE* outputhisto;       // File in which histogram of the density is stored
  double* histogram;       // array in which the density histograms is computed
  int width;               // Histograms is made using squares of width x width boxes
  double rhomax;           // Maximal density used for the histogram
  double drho;             // width of histogram bins
  double NextStoreHisto;   // Next time at which to store the histogram
  double StoreHistoInter;  // Interval between two storage of the histogram
  double HistoInter;       // Interval between two measurement of the histogram
  double NextHisto;        // Next time at which histogram is computed
  int Nbin;                // Number of bins in the density histogram
  double boxarea;          // Area of the box over which densities are computed
  long n0,part;            // counter of particles
  double histocount;       // Number of times a density is recorded
  
  /* SPECIFIC TO RECORDING OF POSITIONS */
  
  double EquilibTimePos; //time after which recordings start for positions
  double StoreInterPos; // Interval at which positions are stored
  double NextStorePos; // next time at which positions are stored
  FILE* outputpos, *outputparam; // Files in which data and parameters are stored
  double* density; //density[i] is the density of particle around particle i
  
#ifdef PRESSURE
  /*
    pressure[0] is the force exerted on the left wall
    pressure[1] is the force exerted on the right wall
  */
  double* pressure = (double*) calloc(sizeof(double),2);
  FILE*   outputpressure;      // File in which histogram of the density is stored
  double  NextStorePressure;   // Next time at which to store the histogram
  double  StorePressureInter;  // Interval between two storage of the histogram
  double  average_density;     // Compute the average density in [Param.Lx*2/5,Param.Lx*3/5]
#endif

#ifdef STRESSTENSOR
  double *** sigma_IK; // array in which the Irving-Kirkwood stress tensor is stored
  double ** sigma_a;  // array in which the active stress tensor is stored
  
  FILE*   output_stress_tensor;       // File in which histogram of the stress tensor is stored
  double  Next_Measure_Stress_Tensor; // Next time at which the stress tensor is measured
  double  Inter_Measure_Stress_Tensor;  // Interval between two storage of the histogram
  double  Next_Store_Stress_Tensor;    // Next time at which the stress tensor is stored
  double  Inter_Store_Stress_Tensor;  // Interval between two storage of the histogram
#endif

#ifdef  GIVENIC
  FILE* input;
#endif

#ifdef SLABIC
  double rhog; // density of particle in the gas phase
  double rhol; // density of particle in the liquid phase
  long Ngas;    // Number of particles in the gas
  long Nliquid; // Number of particles in the liquid
  double liquidfraction; // fraction of system starting in the liquid phase
#endif

#ifdef OBSERVABLE
  double Energy;
  double NextStoreEnergy;
  double StoreInterEnergy;
  FILE* outputenergy;      // File in which histogram of the density is stored
#endif
  
  /* Read parameters of the simulation*/
  
  // Take input from command line, check that their number is correct and store them
  CheckAndReadInput(command_base,argc,argv,&rho0,&Param,&seed,&FinalTime,&EquilibTimePos,&StoreInterPos,&width,&rhomax,&HistoInter,&StoreHistoInter,&outputparam,&outputpos,&outputhisto
#ifdef PRESSURE
		    ,&outputpressure,&StorePressureInter,&average_density
#endif
#ifdef STRESSTENSOR
		    ,&output_stress_tensor,&Inter_Measure_Stress_Tensor,&Inter_Store_Stress_Tensor
#endif
#ifdef GIVENIC
		    ,&input
#endif
#ifdef SLABIC
		    ,&rhog,&rhol,&liquidfraction
#endif
#ifdef OBSERVABLE
		    ,&outputenergy,&StoreInterEnergy
#endif
		    );
  
  /* Initialize variables */
#ifdef _MT
  init_genrand64(seed);
#endif
#ifdef _PCG
  pcg64_srandom(seed,123); 
#endif
  _time               = 0;
  Param.sqrt2Drdt     = sqrt(2*Param.Dr*Param.dt);
  Param.sqrt2Dtdt     = sqrt(2*Param.Dt*Param.dt);
  NextStorePos        = EquilibTimePos;
  
#if defined LJ
  Param.rmax2         = 7.29*Param.sigma*Param.sigma;   // The cutoff is at 2.7*sigma (2.7*2.7=7.29)
#endif

#ifdef WCA
  Param.rmax2         = pow(2,1./3.)*Param.sigma*Param.sigma;   // The cutoff is at 2^{1/6}*sigma
#endif

  Param.sigma2        = pow(Param.sigma,2);

#if defined(WCA) || defined (LJ)
  Param.sigma6        = pow(Param.sigma,6);
  Param.amp           = 24*Param.epsilon*Param.sigma6;
#endif
  
#ifdef RANDOMIC
#ifdef CLOSEDBC
  Param.N             = (long) ( rho0 * (Param.x_wall_right-Param.x_wall_left) *Param.Ly );
#else
  Param.N             = (long) (rho0*Param.Lx*Param.Ly);
#endif
#endif
  
#ifdef SLABIC
  Ngas                = (long) ( rhog * Param.Lx *Param.Ly * (1-liquidfraction) );
  Nliquid             = (long) ( rhol * Param.Lx *Param.Ly * liquidfraction     );
  Param.N             = Ngas+Nliquid;
#endif
  
  Particles           = (particle*) malloc(Param.N*sizeof(particle));
  density             = (double*) calloc(Param.N,sizeof(double));
  forces              = (double*) calloc(Param.N*2,sizeof(double));
  
  // Create boxes in which to store particles for the spatial hashing
  Param.rbox2         = Param.rbox * Param.rbox;
  
  //Parameters for the histogram
  boxarea         = width*width*Param.rbox*Param.rbox;
  Nbin            = (int) (1+floor(rhomax*boxarea));
  histogram       = (double*) calloc(Nbin,sizeof(double));
  NextStoreHisto  = StoreHistoInter;
  NextHisto       = HistoInter;
  histocount      = 0;
  
#ifdef PRESSURE
  NextStorePressure = StorePressureInter;
#endif

#ifdef STRESSTENSOR
  sigma_IK = (double***) malloc(Param.NxBox*sizeof(double*));
  sigma_a  = (double**) malloc(Param.NxBox*sizeof(double*));
  for (i=0;i<Param.NxBox;i++){
    sigma_IK[i] = (double**) calloc(Param.NyBox,sizeof(double));
    sigma_a[i]  = (double*) calloc(Param.NyBox,sizeof(double));
    for (j=0;j<Param.NyBox;j++)
      sigma_IK{i][j] = (double*) calloc(4,sizeof(double));
  }
  
  Next_Measure_Stress_Tensor = Inter_Measure_Stress_Tensor;
  Next_Store_Stress_Tensor   = Inter_Store_Stress_Tensor;
  count_measure_stress_tensor= 0;
#endif

#ifdef OBSERVABLE
  Energy=0;
  NextStoreEnergy   = StoreInterEnergy;
#endif
  
  
#if defined(WCA) || defined(LJ)
  if(Param.rbox2<Param.rmax2)
    ERROR("Box size smaller than interaction length\n");
#endif
  
  Param.NxBox = (long) (floor(Param.Lx/Param.rbox)+EPS);
  if(fabs((double)Param.NxBox*Param.rbox-Param.Lx)>EPS)
    ERROR("Lx has to be a multiple of rbox");
  
  Param.NyBox = (long) (floor(Param.Ly/Param.rbox)+EPS);
  if(fabs((double)Param.NyBox*Param.rbox-Param.Ly)>EPS)
    ERROR("Ly has to be a multiple of rbox");
  
  Box       = (long**) malloc(Param.NxBox*sizeof(long*));
  // Boxes are initialized as empty
  for (i=0;i<Param.NxBox;i++){
    Box[i] = (long*) calloc(Param.NyBox,sizeof(long));
    for (j=0;j<Param.NyBox;j++)
      Box[i][j]=-1;
  }
  
  // This list the neighbours of particles. Initialized as empty.
  Neighbours = (long*) calloc((Param.N+1)*2,sizeof(long));
  for(i=0;i<Param.N;i++){
    Neighbours[2*i]=-1;
    Neighbours[2*i+1]=-1;
  }
  
  //This is a list of all the boxes around box i in which particles
  //may interact with those in box i epsilonx and epsilony are offset
  //that are used for periodic boundary conditions
  NeighbouringBoxes = (box***) malloc(Param.NxBox*sizeof(box**));
  DefineNeighbouringBoxes(NeighbouringBoxes,Param);
  
  /* Setup initial conditions */
  printf("Particle initialisation\n");
  
  // On a lattice
  //Particles[i].x=(i/d)*Param.Lx/d+Param.rmin*genrand64_real3();
  //Particles[i].y=(i%d)*Param.Ly/d+Param.rmin*genrand64_real3();    
  
#ifdef RANDOMIC
  Place_Random_IC_Continuous(Param,Particles,Box,Neighbours);
#endif    
  
#ifdef GIVENIC
  Place_Given_IC(Param,input,Particles,Box,Neighbours);
#endif    
    
#ifdef SLABIC
  Place_Slab_IC_Lattice(Param,liquidfraction,outputparam,Particles,Box,Neighbours,Ngas,Nliquid);
#endif

  // Store Initial conditions
  Density_TopHat_Hashing_PBC(Particles,Param,density,Box,Neighbours,NeighbouringBoxes);
  for(i=0;i<Param.N;i++)
    fprintf(outputpos,"%lg\t%d\t%lg\t%lg\t%lg\t%lg\n",_time,i,Particles[i].x,Particles[i].y,Particles[i].theta,density[i]*Param.sigma2);
  fprintf(outputpos,"\n");
  
  //Store Parameters
  StoreParam(Param,outputparam,command_base,argc, argv,rho0, seed,FinalTime,EquilibTimePos,StoreInterPos,width, rhomax, HistoInter,Nbin,StoreHistoInter
#ifdef PRESSURE
	     , StorePressureInter
#endif
#ifdef STRESSTENSOR
	     , Inter_Measure_Stress_Tensor, Inter_Store_Stress_Tensor
#endif
#ifdef OBSERVABLE
	     , StoreInterEnergy
#endif
	     );
  
  printf("Initialisation over\n");
  
  time_clock       = time(NULL);
  
  /* Run the dynamics */
  
  while(_time<FinalTime){
    // Move the particles
    Move_Particles_ABP(Particles,Param,forces,_time,Box,Neighbours,NeighbouringBoxes,outputparam
#ifdef PRESSURE
		       ,pressure,&average_density
#endif
#ifdef STRESSTENSOR
		       ,sigma_IK,sigma_a,&next_measure_stress_tensor,&next_store_stress_tensor,inter_store_stress_tensor, &count_measure_stress_tensor
#endif
#ifdef OBSERVABLE
		       ,&Energy
#endif
		       );
    
    // Increment time
    _time += Param.dt;
    
    //If it is time, store positions
    if(_time>NextStorePos-EPS)
      StorePositions(Particles, Param, density,Box,Neighbours, NeighbouringBoxes,outputpos,_time,&NextStorePos,StoreInterPos);
    
    //If it is time to update the histogram
    if(_time>NextHisto-EPS)
      UpdateHistogram(Param, width,Neighbours,Nbin,histogram,&histocount,&NextHisto,HistoInter,Box,outputparam,rhomax,_time);

    //if it is time to store the histogram
    if(_time>NextStoreHisto-EPS)
      StoreHistogram(Nbin,outputhisto,_time,boxarea,histogram,&histocount, &NextStoreHisto, StoreHistoInter,Param);
    
#ifdef PRESSURE
    //if it is time to record the pressure
    if(_time>NextStorePressure-EPS)
      StorePressure(outputpressure,  _time, pressure, Param,&average_density,StorePressureInter,&NextStorePressure);
#endif

#ifdef OBSERVABLE
    if(_time>NextStoreEnergy-EPS)
      StoreEnergy(outputenergy,  _time, &Energy, Param,StoreInterEnergy,&NextStoreEnergy);
#endif
  }
  
  printf("Simulation time: %ld seconds\n",time(NULL)-time_clock);
  fprintf(outputparam,"#Simulation time: %ld seconds\n",time(NULL)-time_clock);
  
  free(histogram);
  free(Particles);
  free(density);
  free(forces);
  for (i=0;i<Param.NxBox;i++){
    free(Box[i]);
  }
  free(Box);
  free(Neighbours); 
  for (i=0;i<Param.NxBox;i++){
    for (j=0;j<Param.NyBox;j++)
      free(NeighbouringBoxes[i][j]);
    free(NeighbouringBoxes[i]);
  }
  free(NeighbouringBoxes);
  fclose(outputparam);
  fclose(outputpos);
  fclose(outputhisto);

#ifdef PRESSURE
  fclose(outputpressure);
  free(pressure);
#endif

#ifdef STRESSTENSOR
  for (i=0;i<Param.NxBox;i++){
    for (j=0;j<Param.NxBox;j++)
      free(sigma_IK[i][j]);
    free(sigma_IK[i]);
    free(sigma_a[i]);
  }
  free(sigma_IK);
  free(sigma_a);
  
#endif
  
#ifdef GIVENIC
  fclose(input);
#endif


  return 1;
}


