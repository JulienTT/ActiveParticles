// TODO: boxes and BC in the closed case
// Check what happens when particles exit boxes

// Initial condition: add a slab IC using triangular lattice


/*
  2021-02-24 
  This code simulates the dynamics of active Brownian particles whose dynamics is given by
  
  \dot x = v_0 \cos\theta - \partial_x F + \sqrt{2 D_t} \eta_x 
  \dot y = v_0 \sin\theta - \partial_y F + \sqrt{2 D_t} \eta_y
  \dot \theta = \sqrt{2 D_r} \eta_\theta

  We use Lennard-Jones interactions for F = -\grad V where

  V=epsilon (sigma^12/r^12-sigma^6/r^6)

  We use a cut-off at r=sigma*2.7

  There is for now no adaptative time-stepping because of the order dt and \sqrt{dt}.
  
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mt19937-64.c"
#include <time.h>

/*
  Choose the initial condition between random, if RANDOMIC is defined
  and read from an input if GIVENIC is defined
 */

//#define RANDOMIC
#define GIVENIC


/*
  Choose the boundary conditions between periodic, if PBC is defined,
  and closed along x, if CLOSEDBC is defined
*/
//#define PBC
#define CLOSEDBC

/*
  Uncomment to compute the pressure exerted on right and left walls
*/
#define PRESSURE

#define EPS 1e-10

// Structure particles which contains all the data needed to characterize the state of a particle
typedef struct particle{
  double x;
  double y; //position of particles
  double theta; // orientation of particles
  long bi; // position of box along x axis
  long bj; // position of box along y axis
} particle;

typedef struct box{
  long i; // position of box along x axis
  long j; // position of box along y axis
  double epsilonx; //distance to add along x to take into account PBC
  double epsilony; //distance to add along y to take into account PBC
} box;

// Structure param contains the parameters of the code that will be passed to functions in a compact way
struct param {
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
  double sigma6; //sigma^6
  double rmax2;  // interaction cut-off (squared)
  double amp; //used to compute the force
  
  // Parameters of the spatial hashing
  double rbox; //Box size used for hashing. Equal rmax, a priori.
  double rbox2;
  long NxBox; // number of boxes along x
  long NyBox; // number of boxes along y

#ifdef CLOSEDBC
  /*
    The potential is Omega/nu (x \pm L_x)^nu
  */
  double Omega;
  double nu;
#endif
  
};

#include "ABP-LJ-functions.c"

int main(int argc, char* argv[]){
  
  /* Variable declaration */
  long i,j,k,i1,j1; // counter
  time_t time_clock; // current physical time
  double _time; // current time
  double rho0;
  double nlocal; // number of particle in a box
  
  /*
    To construct a histogram of the data, we average over a period
    StoreHistoInter in an array histogram. We bin together several
    neighbouring boxes, over a linear size width, so that the
    resolution is 1/width^2.

    At the end, we plot rho*sigma^2 so that the resolution is  (sigma/width)^2
  */
  FILE* outputhisto;       // File in which histogram of the density is stored
  double* histogram;       // array in which the density histograms is computed
  double width;            // Histograms is made using squares of width x width boxes
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
  
  double EquilibTimePos; //time after which recordings start for positions
  double StoreInterPos; // Interval at which positions are stored
  double NextStorePos; // next time at which positions are stored
  
  double FinalTime; //final time of the simulation
  
  particle* Particles; // arrays containing the particles
  struct param Param; // Structure containing the parameters
  long long seed; // see of the random number generator
  FILE* outputpos, *outputparam; // Files in which data and parameters are stored
  char name[200]; // string in which the file names are written
  double* density;
  double* forces ; // This array is used to compute the forces exerted
	           // on the particles. forces[2*i] is the force along
	           // x on particle i. forces[2*i+1] is the force
	           // along y on particle i.
  long** Box;            // array containing the first particles in the box
                        // Box[i][j]=k means the first particle in box i,j is k.
  long* Neighbours;   // array containing the neighbours of each
		     // particle. Neighbours[2*i+1]=k means the next
		     // particle after i is k. k=-1 means i is the
		     // last particle. Neighbours[2*i]=k means the
		     // particle before i is k. k=-1 means the particle
		     // is the first in the box.
  box *** NeighbouringBoxes;

#ifdef PRESSURE
  /*
    pressure[0] is the force exerted on the left wall
    pressure[1] is the force exerted on the right wall
  */
  double* pressure = (double*) calloc(sizeof(double),2);
  FILE* outputpressure;       // File in which histogram of the density is stored
  double NextStorePressure;   // Next time at which to store the histogram
  double StorePressureInter;  // Interval between two storage of the histogram
  double PressureInter;       // Interval between two measurement of the histogram
  double NextPressure;        // Next time at which histogram is computed
  double pressurecount;       // Number of times a density is recorded
#endif
  
#ifdef  GIVENIC
  FILE* input;
#endif

  
  /* Read parameters of the simulation*/

  int argctarget=19;
  char command_base[1000]="";
  strcat(command_base, "usage: ");
  strcat(command_base, argv[0]);
  strcat(command_base," file rho0 Lx Ly dt v0 Dr Dt sigma epsilon seed FinalTime EquilibTimePos StoreInterPos width rhomax HistoInter StoreHistoInter");
  
#ifdef GIVENIC
  strcat(command_base," input N");
  argctarget += 2;
#endif
  
#ifdef CLOSEDBC
  strcat(command_base," Omega nu");
  argctarget += 2;
#endif
  
#ifdef PRESSURE
  strcat(command_base," PressureInter StorePressureInter");
  argctarget += 2;
#endif
  if(argc!=argctarget){
    printf("%s\n",command_base);
    exit(1);
  }
  
  i=1;
  //File where parameters are stored
  sprintf(name,"%s-param",argv[i]);
  outputparam=fopen(name,"w");
  
  //File where positions are stored
  sprintf(name,"%s-pos",argv[i]);
  outputpos=fopen(name,"w");

  //File where histogram is stored
  sprintf(name,"%s-histo",argv[i]);
  outputhisto=fopen(name,"w");

#ifdef PRESSURE
  //File where histogram is stored
  sprintf(name,"%s-pressure",argv[i]);
  outputpressure=fopen(name,"w");
#endif
  
  i++;
  rho0             = strtod(argv[i], NULL); i++;  
  Param.Lx         = strtod(argv[i], NULL); i++;
  Param.Ly         = strtod(argv[i], NULL); i++;
  Param.dt         = strtod(argv[i], NULL); i++;
  Param.v0         = strtod(argv[i], NULL); i++;
  Param.Dr         = strtod(argv[i], NULL); i++;
  Param.Dt         = strtod(argv[i], NULL); i++;
  Param.sigma      = strtod(argv[i], NULL); i++;
  Param.epsilon    = strtod(argv[i], NULL); i++; 
  seed             = (long long) strtod(argv[i], NULL); i++;
  FinalTime        = strtod(argv[i], NULL); i++;
  EquilibTimePos   = strtod(argv[i], NULL); i++;
  StoreInterPos    = strtod(argv[i], NULL); i++;
  width            = strtod(argv[i], NULL); i++;
  rhomax           = strtod(argv[i], NULL); i++;
  HistoInter       = strtod(argv[i], NULL); i++;
  StoreHistoInter  = strtod(argv[i], NULL); i++;

#ifdef GIVENIC
  sprintf(name,"%s",argv[i]);i++;
  printf("read from file %s\n",name);
  Param.N          = (long) strtod(argv[i], NULL); i++;
  printf("N is %ld\t rho is %lg\n",Param.N,Param.N/Param.Lx/Param.Ly);
  input            = fopen(name,"r");
#endif

#ifdef CLOSEDBC
  Param.Omega      = strtod(argv[i], NULL); i++;
  Param.nu         = strtod(argv[i], NULL); i++;
#endif
  
#ifdef PRESSURE
  PressureInter      = strtod(argv[i], NULL); i++;
  StorePressureInter = strtod(argv[i], NULL); i++;
#endif
  
  /* Initialize variables */
  init_genrand64(seed);
  _time               = 0;
  
  Param.sqrt2Drdt     = sqrt(2*Param.Dr*Param.dt);
  Param.sqrt2Dtdt     = sqrt(2*Param.Dt*Param.dt);
  
  NextStorePos        = EquilibTimePos;
  
  Param.rmax2         = 7.29*Param.sigma*Param.sigma;   // The cutoff is at 2.7*sigma (2.7*2.7=7.29)
  Param.sigma2        = pow(Param.sigma,2);
  Param.sigma6        = pow(Param.sigma,6);
  Param.amp           = 24*Param.epsilon*Param.sigma6;

#ifdef RANDOMIC
  Param.N             = (long) (rho0*Param.Lx*Param.Ly);
#endif
  
  Particles           = (particle*) malloc(Param.N*sizeof(particle));
  density             = (double*) malloc(Param.N*sizeof(double));
  forces              = (double*) malloc(Param.N*2*sizeof(double));

  // Create boxes in which to store particles for the spatial hashing
  Param.rbox=1;
  Param.rbox2=  Param.rbox*  Param.rbox;

  //Parameters for the histogram
  boxarea         = width*width*Param.rbox*Param.rbox;
  Nbin            = (int) (1+floor(rhomax*boxarea));
  histogram       = (double*) calloc(Nbin,sizeof(double));
  NextStoreHisto  = StoreHistoInter;
  NextHisto       = HistoInter;
  histocount      = 0;

#ifdef PRESSURE
  NextStorePressure = StorePressureInter;
  NextPressure      = PressureInter;
  pressurecount     = 0;
#endif
  
  if(Param.rbox2<Param.rmax2){
    printf("Box size smaller than interaction length\n");
    exit(1);
  }
  
  Param.NxBox = (long) (floor(Param.Lx/Param.rbox)+EPS);
  if(fabs((double)Param.NxBox*Param.rbox-Param.Lx)>EPS){
    printf("Lx has to be a multiple of rbox");
    exit(1);
  }
  
  Param.NyBox = (long) (floor(Param.Ly/Param.rbox)+EPS);
  if(fabs((double)Param.NyBox*Param.rbox-Param.Ly)>EPS){
    printf("Ly has to be a multiple of rbox");
    exit(1);
  }
  
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
  for (i=0;i<Param.NxBox;i++){
    NeighbouringBoxes[i] = (box**) malloc(Param.NyBox*sizeof(box*));
    for (j=0;j<Param.NyBox;j++){
      NeighbouringBoxes[i][j] = (box*) malloc(5*sizeof(box));
      
      // This is the same box
      NeighbouringBoxes[i][j][0].i=i;
      NeighbouringBoxes[i][j][0].j=j;
      NeighbouringBoxes[i][j][0].epsilonx=0;
      NeighbouringBoxes[i][j][0].epsilony=0;
      
      // The box above
      NeighbouringBoxes[i][j][1].i=i;
      NeighbouringBoxes[i][j][1].j=(j+1)%Param.NyBox;
      NeighbouringBoxes[i][j][1].epsilonx=0;
      NeighbouringBoxes[i][j][1].epsilony=(j==Param.NyBox-1)?(Param.Ly):0;

      // The box above, to the right
      NeighbouringBoxes[i][j][2].i=(i+1)%Param.NxBox;
      NeighbouringBoxes[i][j][2].j=(j+1)%Param.NyBox;
      NeighbouringBoxes[i][j][2].epsilonx=(i==Param.NxBox-1)?(Param.Lx):0;
      NeighbouringBoxes[i][j][2].epsilony=(j==Param.NyBox-1)?(Param.Ly):0;

      // The box to the right
      NeighbouringBoxes[i][j][3].i=(i+1)%Param.NxBox;
      NeighbouringBoxes[i][j][3].j=j;
      NeighbouringBoxes[i][j][3].epsilonx=(i==Param.NxBox-1)?(Param.Lx):0;
      NeighbouringBoxes[i][j][3].epsilony=0;

      // The box below, to the right
      NeighbouringBoxes[i][j][4].i=(i+1)%Param.NxBox;
      NeighbouringBoxes[i][j][4].j=(j-1+Param.NyBox)%Param.NyBox;
      NeighbouringBoxes[i][j][4].epsilonx=(i==Param.NxBox-1)?(Param.Lx):0;
      NeighbouringBoxes[i][j][4].epsilony=(j==0)?(-Param.Ly):0;
    }
  }

  /* Setup initial conditions */
  printf("Particle initialisation\n");
  
#ifdef RANDOMIC
  for(i=0;i<Param.N;i++){
    double mindist2,distance2;
    /*I.C. : here the particles are disposed according to a grid. But
      It's really regular only if N is a square !*/
    
    //Particles[i].x=(i/d)*Param.Lx/d+Param.rmin*genrand64_real3();
    //Particles[i].y=(i%d)*Param.Ly/d+Param.rmin*genrand64_real3();    
    
    /* 
       Here we uniformly randomly draw particles one by one, rejecting
       a particle's initial position if it's too close to another.To
       be used for hard core interaction.
    */
    mindist2 = 0;//initialization so that we pass at least once in the loop.
    
    while(mindist2<Param.sigma2){//Critère un peu arbitraire...
      mindist2=Param.Lx*Param.Ly;
      Particles[i].x  = Param.Lx*genrand64_real3();
      Particles[i].y  = Param.Ly*genrand64_real3();
      for(j=0;j<i;j++){
	distance2 = Distance2(Particles,i,j,Param);
	if(distance2<mindist2)
	  mindist2=distance2;
      }
      //printf("particle %ld, min distance %lg\n",i,sqrt(mindist2));
    }
    
    Particles[i].theta = 2*M_PI*genrand64_real2();
    
    //Add particle in the good box.
    Particles[i].bi    = floor(Particles[i].x/Param.rbox);
    Particles[i].bj    = floor(Particles[i].y/Param.rbox);
    AddinBox(i,Particles[i].bi,Particles[i].bj,Box,Neighbours);
  }
  printf("Smallest distance %lg\n",SmallestDistance(Particles,Param));
#endif    
  
#ifdef GIVENIC
  for(i=0;i<Param.N;i++){
    double x,y,theta;
    if(fscanf(input,"%*lg\t%*d\t%lg\t%lg\t%lg\t%*lg\n",&x,&y,&theta)!=EOF){
      Particles[i].x=x;
      Particles[i].y=y;
      Particles[i].theta=theta;
      //printf("Particle %ld at %lg %lg with angle %lg\n",i,x,y,theta);
      
      //Add particle in the good box.
      Particles[i].bi    = floor(Particles[i].x/Param.rbox);
      Particles[i].bj    = floor(Particles[i].y/Param.rbox);
      AddinBox(i,Particles[i].bi,Particles[i].bj,Box,Neighbours);
    }
    else{
      printf("Out of lines\n");
      exit(1);
    }
  }
#endif    
  
  /* Store parameters */
  fprintf(outputparam,"usage: %s file rho0 Lx Ly dt v0 Dr Dt sigma epsilon seed FinalTime EquilibTimePos StoreInterPos width rhomax HistoInter StoreHistoInter\n",argv[0]);
  
  for(i=0;i<argc;i++){
    fprintf(outputparam,"%s ",argv[i]);
  }
  fprintf(outputparam,"\n");
  
  fprintf(outputparam,"FinalTime is %lg\n",FinalTime);
  fprintf(outputparam,"rho0 is %lg\n", rho0);
  fprintf(outputparam,"N is %ld\n", Param.N);
  fprintf(outputparam,"Lx is %lg\n", Param.Lx);
  fprintf(outputparam,"Ly is %lg\n", Param.Ly);
  fprintf(outputparam,"dt is %lg\n", Param.dt);
  fprintf(outputparam,"v0 is %lg\n", Param.v0);
  fprintf(outputparam,"Dt is %lg\n", Param.Dt);
  fprintf(outputparam,"Dr is %lg\n", Param.Dr);
  fprintf(outputparam,"sigma is %lg\n", Param.sigma);
  fprintf(outputparam,"epsilon is %lg\n", Param.epsilon);
  fprintf(outputparam,"seed is %lld\n", seed);
  fprintf(outputparam,"EquilibTimePos is %lg\n", EquilibTimePos);
  fprintf(outputparam,"StoreInterPos is %lg\n", StoreInterPos);
  fprintf(outputparam,"width is %lg\n", width);
  fprintf(outputparam,"rhomax is %lg\n", rhomax);
  fprintf(outputparam,"histointer is %lg\n", HistoInter );
  fprintf(outputparam,"storehistointer is %lg\n", StoreHistoInter);
  fprintf(outputparam,"nbins is %lg\n", Nbin);
  fprintf(outputparam,"rmax2 is %lg\n", Param.rmax2);
  fprintf(outputparam,"rbox2 is %lg\n", Param.rbox2);

  fprintf(outputparam,"normalized packing fraction is N sigma^2 /(Lx Ly) = %lg\n", Param.N*Param.sigma*Param.sigma/Param.Lx/Param.Ly);
  
  fflush(outputparam);
  
  // Save initial condition
  /*  Density_TopHat_Hashing_PBC(Particles,Param,density,Box,Neighbours,NeighbouringBoxes); */
  /*  Kernel_TwoSchwartz_Hashing_PBC(Particles,Param,velocity,Box,Neighbours,NeighbouringBoxes); */
  
  /* for(i=0;i<Param.N;i++){ */
  /*   fprintf(output,"%lg\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\n",_time,i,Particles[i].x,Particles[i].y,Particles[i].theta,density[i],Param.v0*exp(Param.beta*velocity[i])); */
  /* } */
  /* fprintf(output,"\n"); */
  /* fflush(output); */
  
  printf("Initialisation over\n");
  
  time_clock       = time(NULL);
  
  /* Run the dynamics */
  
  while(_time<FinalTime){
    // Move the particles
#ifdef PRESSURE
    Move_Particles_ABP(Particles,Param,forces,_time,Box,Neighbours,NeighbouringBoxes,outputparam,pressure,&pressurecount);
#else
    Move_Particles_ABP(Particles,Param,forces,_time,Box,Neighbours,NeighbouringBoxes,outputparam);
#endif
    
    // Increment time
    _time += Param.dt;
    
    //If it is time, store positions
    if(_time>NextStorePos-EPS){
      //Compute density
      Density_TopHat_Hashing_PBC(Particles,Param,density,Box,Neighbours,NeighbouringBoxes);
      
      // Store positions, angles, packing fraction using a particle diameter sigma
      for(i=0;i<Param.N;i++){
	fprintf(outputpos,"%lg\t%d\t%lg\t%lg\t%lg\t%lg\n",_time,i,Particles[i].x,Particles[i].y,Particles[i].theta,density[i]*Param.sigma2);
	//fprintf(outputpos,"%lg\t%d\t%lg\t%lg\t%lg\n",_time,i,Particles[i].x,Particles[i].y,Particles[i].theta);
      }
      fprintf(outputpos,"\n");
      fflush(outputpos);
      
      //Increase NextStore
      NextStorePos += StoreInterPos;
    }
    
    //If it is time to record the histogram
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
      
      
      //if it is time to record the histogram
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
      fprintf(outputpressure,"%lg\t%lg\t%lg\n",_time,pressure[0]/pressurecount,pressure[1]/pressurecount);
      NextStorePressure += StorePressureInter;
      pressurecount=0;
      pressure[0]=0;
      pressure[1]=0;
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
  
  return 1;
}


