// TODO: boxes and BC in the closed case
// Check what happens when particles exit boxes

// Initial condition: add a slab IC using triangular lattice

/*
  2021-02-24 
  This code simulates the dynamics of active Brownian particles whose dynamics is given by
  
  \dot x = v_0 \cos\theta - \partial_x F + \sqrt{2 D_t} \eta_x 
  \dot y = v_0 \sin\theta - \partial_y F + \sqrt{2 D_t} \eta_y
  \dot \theta = \sqrt{2 D_r} \eta_\theta
*/


//#define LJ
/*
  We use Lennard-Jones interactions for F = -\grad V where
  
  V=4 epsilon (sigma^12/r^12-sigma^6/r^6)
  
  We use a cut-off at r=sigma*2.7
  
  There is for now no adaptative time-stepping because of the order dt and \sqrt{dt}.
  
*/

#define WCA
/*
  We use WCA potential for F = -\grad V where
  
  V=4 epsilon (sigma^12/r^12-sigma^6/r^6)
  
  We use a cut-off at r=sigma*2^1/6
  
  There is for now no adaptative time-stepping because of the order dt and \sqrt{dt}.
  
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mt19937-64.c"
#include <time.h>

/*
  Choose the initial condition (IC) between random IC, if RANDOMIC is
  defined, or read from an input if GIVENIC is defined
*/

#define RANDOMIC
//#define GIVENIC

/*
  Choose the boundary conditions between periodic, if PBC is defined,
  and closed along x, if CLOSEDBC is defined
  
  If closed boundary conditions are used, then a confining potential
  is inserted INSIDE the box, to keep the particles far away from the
  periodic boundary condition of the cell. If particles get close
  enough to the wall, the program protests.
  
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
  double mu; //particle mobility
  
  // Parameters of the spatial hashing
  double rbox; //Box size used for hashing. Equal rmax, a priori.
  double rbox2;
  long NxBox; // number of boxes along x
  long NyBox; // number of boxes along y

#ifdef CLOSEDBC
  /*
    The potential is Omega/nu (x_wall_left-x)^nu and (x-x_wall_right)^nu
  */
  double x_wall_left;
  double x_wall_right;
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

  /* OTHER VARIABLES */
  
  double FinalTime; //final time of the simulation
  particle* Particles; // arrays containing the particles
  struct param Param; // Structure containing the parameters
  long long seed; // see of the random number generator
  char name[200]; // string in which the file names are written

  
  double* forces ;   // This array is used to compute the forces exerted
	             // on the particles. forces[2*i] is the force along
	             // x on particle i. forces[2*i+1] is the force
	             // along y on particle i.
  
  long** Box;        // array containing the first particles in the box
                     // Box[i][j]=k means the first particle in box i,j is k.
  
  long* Neighbours;  // array containing the neighbours of each
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
  double  Energypot=0;         // Use to compute the average potential energy
  FILE*   outputpressure;      // File in which histogram of the density is stored
  double  NextStorePressure;   // Next time at which to store the histogram
  double  StorePressureInter;  // Interval between two storage of the histogram
  double  average_density;     // Compute the average density in [Param.Lx*2/5,Param.Lx*3/5]
#endif
  
#ifdef  GIVENIC
  FILE* input;
#endif
  
  /* Read parameters of the simulation*/
  
  int argctarget=21;
  char command_base[1000]="";
  strcat(command_base, "usage: ");
  strcat(command_base, argv[0]);
  strcat(command_base," file rho0 Lx Ly dt v0 Dr Dt mu sigma epsilon seed FinalTime EquilibTimePos StoreInterPos rbox width rhomax HistoInter StoreHistoInter");
  
#ifdef GIVENIC
  strcat(command_base," input N");
  argctarget += 2;
#endif
  
#ifdef CLOSEDBC
  strcat(command_base," Omega nu x_wall_left x_wall_right");
  argctarget += 4;
#endif
  
#ifdef PRESSURE
  strcat(command_base," StorePressureInter");
  argctarget += 1;
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
  Param.mu         = strtod(argv[i], NULL); i++;
  Param.sigma      = strtod(argv[i], NULL); i++;
  Param.epsilon    = strtod(argv[i], NULL); i++; 
  seed             = (long long) strtod(argv[i], NULL); i++;
  FinalTime        = strtod(argv[i], NULL); i++;
  EquilibTimePos   = strtod(argv[i], NULL); i++;
  StoreInterPos    = strtod(argv[i], NULL); i++;
  Param.rbox       = strtod(argv[i], NULL); i++;
  width            = (int) strtod(argv[i], NULL); i++;
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
  Param.Omega       = strtod(argv[i], NULL); i++;
  Param.nu          = strtod(argv[i], NULL); i++;
  Param.x_wall_left = strtod(argv[i], NULL); i++;
  Param.x_wall_right= strtod(argv[i], NULL); i++;
#endif
  
#ifdef PRESSURE
  StorePressureInter = strtod(argv[i], NULL); i++;
  average_density    = 0;
#endif
  
  /* Initialize variables */
  init_genrand64(seed);
  _time               = 0;
  
  Param.sqrt2Drdt     = sqrt(2*Param.Dr*Param.dt);
  Param.sqrt2Dtdt     = sqrt(2*Param.Dt*Param.dt);
  
  NextStorePos        = EquilibTimePos;

#ifdef LJ
  Param.rmax2         = 7.29*Param.sigma*Param.sigma;   // The cutoff is at 2.7*sigma (2.7*2.7=7.29)
#endif
#ifdef WCA
  Param.rmax2         = pow(2,1./3.)*Param.sigma*Param.sigma;   // The cutoff is at 2.7*sigma (2.7*2.7=7.29)
#endif
  Param.sigma2        = pow(Param.sigma,2);
  Param.sigma6        = pow(Param.sigma,6);
  Param.amp           = 24*Param.epsilon*Param.sigma6;
  
#ifdef RANDOMIC
#ifdef CLOSEDBC
  Param.N             = (long) ( rho0 * (Param.x_wall_right-Param.x_wall_left) *Param.Ly );
#else
  Param.N             = (long) (rho0*Param.Lx*Param.Ly);
#endif
#endif
  
  Particles           = (particle*) malloc(Param.N*sizeof(particle));
  density             = (double*) calloc(Param.N,sizeof(double));
  forces              = (double*) calloc(Param.N*2,sizeof(double));
  
  // Create boxes in which to store particles for the spatial hashing
  Param.rbox2 = Param.rbox * Param.rbox;
  
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
    printf("Ly has to be a multiple of rbox, right now %lg / %lg = %lg",Param.Ly,Param.rbox,Param.Ly/Param.rbox);
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

  // On a lattice
  //Particles[i].x=(i/d)*Param.Lx/d+Param.rmin*genrand64_real3();
  //Particles[i].y=(i%d)*Param.Ly/d+Param.rmin*genrand64_real3();    

#ifdef RANDOMIC
  /*
    This should be replaced by a rejection method on a triangular
    lattice of length .9 sigma.
   */
  for(i=0;i<Param.N;i++){
    double mindist2,distance2;
    
    /* 
       Here we uniformly randomly draw particles one by one, rejecting
       a particle's initial position if it's too close to another.To
       be used for hard core interaction.
    */
    mindist2 = 0;//initialization so that we pass at least once in the loop.
    
    while(mindist2<.8*Param.sigma2){//Critère un peu arbitraire...
      mindist2=Param.Lx*Param.Ly;

#ifdef CLOSEDBC
      Particles[i].x = Param.x_wall_left+(Param.x_wall_right-Param.x_wall_left)*genrand64_real3();
#else
      Particles[i].x = Param.Lx*genrand64_real3();
#endif
      Particles[i].y = Param.Ly*genrand64_real3();
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
      printf("Not the right number of lines\n");
      exit(1);
    }
  }
#endif    
  
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
      average_density=0;
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

#ifdef PRESSURE
  fclose(outputpressure);
  free(pressure);
#endif

#ifdef GIVENIC
  fclose(input);
#endif
    
  return 1;
}


