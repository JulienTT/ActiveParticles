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
#include "mt19937-64.c"
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
#if defined(WCA) || defined(LJ)
  double epsilon; // amplitude of the potential
  double sigma; //translational diffusivity
  double sigma2; //sigma^2
  double sigma6; //sigma^6
  double rmax2;  // interaction cut-off (squared)
  double amp; //used to compute the force
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

#include "ABP-LJ-functions.c"

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
  long long seed;      // see of the random number generator

  
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
  double  Energypot=0;         // Use to compute the average potential energy
  FILE*   outputpressure;      // File in which histogram of the density is stored
  double  NextStorePressure;   // Next time at which to store the histogram
  double  StorePressureInter;  // Interval between two storage of the histogram
  double  average_density;     // Compute the average density in [Param.Lx*2/5,Param.Lx*3/5]
#endif
  
#ifdef  GIVENIC
  FILE* input;
#endif

#ifdef SLABIC
  double rhog; // density of particle in the gas phase
  double rhol; // density of particle in the liquid phase
  double liquidfraction; // fraction of system starting in the liquid phase
#endif
  
  /* Read parameters of the simulation*/

  // Take input from command line, check that their number is correct and store them
  CheckAndReadInput(command_base,argc,argv,&rho0,&Param,&seed,&FinalTime,&EquilibTimePos,&StoreInterPos,&width,&rhomax,&HistoInter,&StoreHistoInter,&outputparam,&outputpos,&outputhisto
#ifdef PRESSURE
		    ,&outputpressure,&StorePressureInter,&average_density
#endif
#ifdef GIVENIC
		    ,&input
#endif
#ifdef SLABIC
		    ,&rhog,&rhol,&liquidfraction
#endif
		    );
  
  /* Initialize variables */
  init_genrand64(seed);
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
#if defined(WCA) || defined (LJ)
  Param.sigma2        = pow(Param.sigma,2);
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
  
  if(Param.rbox2<Param.rmax2)
    ERROR("Box size smaller than interaction length\n");
  
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
  Place_Slab_IC_Lattice();
#endif
  
  for(i=0;i<Param.N;i++){
    fprintf(outputpos,"%lg\t%d\t%lg\t%lg\t%lg\n",_time,i,Particles[i].x,Particles[i].y,Particles[i].theta);
  }
  fprintf(outputpos,"\n");
  
  /* Store parameters */
  fprintf(outputparam,"%s\n",command_base);
  
  for(i=0;i<argc;i++){
    fprintf(outputparam,"%s ",argv[i]);
  }
  fprintf(outputparam,"\n");
  
  fprintf(outputparam,"rho0 is %lg\n", rho0);
  fprintf(outputparam,"Lx is %lg\n", Param.Lx);
  fprintf(outputparam,"Ly is %lg\n", Param.Ly);
  fprintf(outputparam,"dt is %lg\n", Param.dt);
  fprintf(outputparam,"v0 is %lg\n", Param.v0);
  fprintf(outputparam,"Dr is %lg\n", Param.Dr);
  fprintf(outputparam,"Dt is %lg\n", Param.Dt);
  fprintf(outputparam,"sigma is %lg\n", Param.sigma);
  fprintf(outputparam,"epsilon is %lg\n", Param.epsilon);
  fprintf(outputparam,"seed is %lld\n", seed);
  fprintf(outputparam,"FinalTime is %lg\n",FinalTime);
  fprintf(outputparam,"EquilibTimePos is %lg\n", EquilibTimePos);
  fprintf(outputparam,"StoreInterPos is %lg\n", StoreInterPos);
  fprintf(outputparam,"width is %d\n", width);
  fprintf(outputparam,"rhomax is %lg\n", rhomax);
  fprintf(outputparam,"histointer is %lg\n", HistoInter );
  fprintf(outputparam,"storehistointer is %lg\n", StoreHistoInter);
  fprintf(outputparam,"rbox is %lg\n", Param.rbox);
  
#ifdef CLOSEDBC
  fprintf(outputparam,"Omega is %lg\n", Param.Omega);
  fprintf(outputparam,"nu is %lg\n", Param.nu);
  fprintf(outputparam,"x_wall_left is %lg\n", Param.x_wall_left);
  fprintf(outputparam,"x_wall_right is %lg\n", Param.x_wall_right);
#endif

#ifdef PRESSURE
  fprintf(outputparam,"StorePressureInter is %lg\n",StorePressureInter);
#endif
  
  fprintf(outputparam,"nbins is %lg\n", Nbin);
  fprintf(outputparam,"rmax2 is %lg\n", Param.rmax2);
  fprintf(outputparam,"rbox2 is %lg\n", Param.rbox2);
  fprintf(outputparam,"N is %ld\n", Param.N);
  fprintf(outputparam,"normalized packing fraction is N sigma^2 /(Lx Ly) = %lg\n", Param.N*Param.sigma*Param.sigma/Param.Lx/Param.Ly);
  fflush(outputparam);
  
  printf("Initialisation over\n");
  
  time_clock       = time(NULL);
  
  /* Run the dynamics */
  
  while(_time<FinalTime){
    // Move the particles
#ifdef PRESSURE
    Move_Particles_ABP(Particles,Param,forces,_time,Box,Neighbours,NeighbouringBoxes,outputparam,pressure,&average_density);
#else
    Move_Particles_ABP(Particles,Param,forces,_time,Box,Neighbours,NeighbouringBoxes,outputparam);
#endif
    
    // Increment time
    _time += Param.dt;
    
    //If it is time, store positions
    if(_time>NextStorePos-EPS){
      //Compute density
      Density_TopHat_Hashing_PBC(Particles,Param,density,Box,Neighbours,NeighbouringBoxes);
      
      // Store positions, angles, packing-fraction using a particle diameter sigma
      for(i=0 ; i<Param.N ; i++){
	fprintf(outputpos,"%lg\t%d\t%lg\t%lg\t%lg\t%lg\n",_time,i,Particles[i].x,Particles[i].y,Particles[i].theta,density[i]*Param.sigma2);
	//fprintf(outputpos,"%lg\t%d\t%lg\t%lg\t%lg\n",_time,i,Particles[i].x,Particles[i].y,Particles[i].theta);
      }
      fprintf(outputpos,"\n");
      fflush(outputpos);
      
      //Increase NextStore
      NextStorePos += StoreInterPos;
    }
    
    //If it is time to update the histogram
    if(_time>NextHisto-EPS){
      // Loop through all area of size width
      for(i=0;i<Param.NxBox;i+=width){
	for(j=0;j<Param.NyBox;j+=width){
	  // initialise the number or particles in this box to zero
	  n0=0;
	  // loop through all the boxes in the area
	  for(i1=i;i1<i+width;i1++){
	    for(j1=j;j1<j+width;j1++){
	      // Find the first particle in the box
	      part=Box[i1][j1];
	      // As long as there is a next particle, increment n0
	      while(part != -1){
		n0++;
		part=Neighbours[2*part+1];
	      }
	    }
	  }
	  // n0 is now the total number of particle in an area of size width*r0 * width*r0;
	  if(n0<Nbin){
	    histogram[n0] += 1.;
	    histocount += 1;
	  }
	  else{
	    printf("rhomax too small, n0 is %d and Nbin is %d\n",n0,Nbin);
	    fprintf(outputparam,"rhomax= %lg too small at time %lg\n",rhomax,_time);
	    exit(1);
	  }
	}
      }
      
      //if it is time to store the histogram
      if(_time>NextStoreHisto-EPS){
	for(i=0;i<Nbin;i++){
	  fprintf(outputhisto,"%lg\t%lg\t%lg\n",_time,(double) i/boxarea*Param.sigma2,histogram[i]/histocount*boxarea/Param.sigma2);
	}
	NextStoreHisto += StoreHistoInter;
	memset(histogram, 0, sizeof(double)*Nbin);
	histocount=0;
	fprintf(outputhisto,"\n");
	fflush(outputhisto);
      }
      
      NextHisto += HistoInter;
    }
    
#ifdef PRESSURE
    //if it is time to record the pressure
    if(_time>NextStorePressure-EPS){
      // Save the bulk density in param
      fprintf(outputpressure,"%lg\t%lg\t%lg\t%lg\n",_time,pressure[0]/StorePressureInter/Param.Ly,pressure[1]/StorePressureInter/Param.Ly,average_density/Param.Lx*5/Param.Ly/StorePressureInter);
      
      NextStorePressure += StorePressureInter;
      average_density    = 0;
      pressure[0]        = 0;
      pressure[1]        = 0;
      fflush(outputpressure);
    }
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

#ifdef GIVENIC
  fclose(input);
#endif


  return 1;
}


