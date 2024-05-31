/*
  2021-03-09 
  This code simulates the dynamics of N run-and-tumble particles whose dynamics is given by
  
  \dot x = v_0 \cos\theta - \partial_x F
  \dot y = v_0 \sin\theta - \partial_y F
  
  \theta -> theta' at rate alpha

  F represents obstacles. In practicte, we may not compute the force
  exerted by the obstacles but instead cancel the normal component of
  the velocity.
  
  TO BE IMPLEMENTED.

  Another move function which represents obstacles as potentials in
  case we want overlapping or dense obstacles.

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mt19937-64.c"
#include <time.h>

#include <fftw3.h>

#define PERIODICBC

// Structure particles which contains all the data needed to characterize the state of a particle
typedef struct particle{
  double x; //position of particles along x axis
  double y; //position of particles along y axis
  double theta; // orientation of particles
  double next_tumble; //time at which the next tumble occurs
  long bi; // position of box along x axis
  long bj; // position of box along y axis
} particle;

// Structure param contains the parameters of the code that will be passed to functions in a compact way
typedef struct param{
  long N;  // number of particles
  double Lx; // system size
  double Ly;
  double dt; // time-step
  double v0; // particle speed
  double alpha; // tumbling rate
  long Nobs; // number of obstacles
  double sigma; // radius of the obstacles
  double sigma2; // radius of the obstacles
  
  // Parameters of the spatial hashing
  double rbox; //Box size used for hashing. Equal rmax, a priori.
  double rbox2;
  long NxBox; // number of boxes along x
  long NyBox; // number of boxes along y
  int NxBoxDensity; // number of boxes along x to record density
  int NyBoxDensity; // number of boxes along x to record density
  double dr; // bin width to record density
} param;

// This structure is used to store the boxes that have to been looped
// through by particles. It contains both the index of the box (i,j)
// and the offset which have to be added to the obstacles in this box
typedef struct box{
  long i; // position of box along x axis
  long j; // position of box along y axis
  double epsilonx; //distance to add along x to take into account PBC
  double epsilony; //distance to add along y to take into account PBC
} box;

typedef struct obstacle{
  double x; // position of obstacle along x axis
  double y; // position of obstacle along y axis
  long bi;  // position of box along x axis
  long bj;  // position of box along y axis
  int prev; // previous obstacle in the box
  int next; // next obstacle in the box
} obstacle;
  
#include "RTP-Obstacles-functions.c"
#define EPS 1e-10
  
int main(int argc,char* argv[]){
  
  /* DECLARATION OF VARIABLES */
  particle* Particles;
  long i,j;             // counter
  long long seed;       // seed of the random number generator
  double _time;         // current time of the simulation
  double totaltime;     // total time of the simulation
  param Param;          // A bunch of parameters which are needed by the functions


  int** Box;           // array containing the first obstacle in the box
                        // Box[i][j]=k means the first obstacle in box i,j is k.

  box *** NeighbouringBoxes; // NeighbouringBoxes[i][j][k] is one of
			     // the 9 boxes that have to be looped
			     // through to know whether particles in
			     // box i,j interact with an obstacle
  
  obstacle* Obstacles;  // array containing the obstacles

  double** density;               //density[i][j] is density of particles in [i.dr,(i+1).dr] x [j.dr,(j+1).dr]
  double*** current;
  /*current[i][j][0/1] contains the x/y components of the particle current at r=(i*dr,j*dr)*/
  double densitycount;            // counts the number of time density and currents have been recorded
  
  double InterDensity;            // density & current are updated every InterDensity
  double NextDensity      ;       // Next time at which density & current are updated
  double StoreInterDensity;       // Interval between two recordings of the density & current in a file
  double NextStoreDensity;        // Next time at which density and currents are stored
  double EquilibTimeDensity;      // Interval during which no recording or measurement of density is made
  FILE* outputdens;               // File where currents and densities are stored
  
  FILE* outputpos;   // where positions of particles are stored
  FILE* outputparam; // where positions of particles are stored
  FILE* outputobs;   // where parameters of simulations are stored
  char name[100];    // String where file names is written
  
  double EquilibTimePos; //time after which recordings start for positions
  double StoreInterPos; // Interval at which positions are stored
  double NextStorePos; // next time at which positions are stored

  /* INITIALISATION OF VARIABLES */

  int argctarget=18;
  char command_base[1000]="";
  strcat(command_base, "usage: ");
  strcat(command_base, argv[0]);
  strcat(command_base," outputfiles N seed totaltime Lx Ly alpha v0 dt Nobstacles sigma EquilibTimePos StoreInterPos EquilibTimeDensity InterDensity StoreInterDensity dr");
  
  if(argc!=argctarget){
    printf("%s\n",command_base);
    exit(1);
  }
  // Create an array of N particles
  i=1;
  
  //File where positions of particles are stored
  sprintf(name,"%s-pos",argv[i]);
  outputpos=fopen(name,"w");

  //File where positions of obstacles are stored
  sprintf(name,"%s-obs",argv[i]);
  outputobs=fopen(name,"w");

  //File where density and current are stored
  sprintf(name,"%s-dens",argv[i]);
  outputdens=fopen(name,"w");

  //File where parameters are stored
  sprintf(name,"%s-param",argv[i]);
  outputparam=fopen(name,"w");

  
  i++;
  
  Param.N           = (long) strtod(argv[i], NULL); i++;
  seed              = (long long) strtod(argv[i], NULL); i++;
  totaltime         = strtod(argv[i], NULL); i++;
  Param.Lx          = strtod(argv[i], NULL); i++;
  Param.Ly          = strtod(argv[i], NULL); i++;
  Param.alpha       = strtod(argv[i], NULL); i++;
  Param.v0          = strtod(argv[i], NULL); i++;
  Param.dt          = strtod(argv[i], NULL); i++;
  Param.Nobs        = (long) strtod(argv[i], NULL); i++;
  Param.sigma       = strtod(argv[i], NULL); i++;
  EquilibTimePos    = strtod(argv[i], NULL); i++;
  StoreInterPos     = strtod(argv[i], NULL); i++;
  EquilibTimeDensity=strtod(argv[i], NULL); i++;
  InterDensity      =strtod(argv[i], NULL); i++;
  StoreInterDensity =strtod(argv[i], NULL); i++;
  Param.dr          = strtod(argv[i], NULL); i++;
  
  Param.sigma2      = Param.sigma*Param.sigma;
  Param.rbox        = 1;
  Param.rbox2       = 1;

  if(Param.rbox<Param.sigma){
    printf("The obstacles are too big for the grid\n");
    exit(1);
  }
  
  init_genrand64(seed);
  Particles         = (particle*) malloc(sizeof(particle) * Param.N);
  Obstacles         = (obstacle*) malloc(sizeof(obstacle) * Param.Nobs);
  _time             = 0;
  NextStorePos      = EquilibTimePos;
  NextDensity       = EquilibTimeDensity;
  NextStoreDensity  = EquilibTimeDensity+StoreInterDensity;
  
  Param.NxBoxDensity      = Param.Lx/Param.dr+1;
  Param.NyBoxDensity      = Param.Ly/Param.dr+1;
  
  density           = (double**) malloc(sizeof(double*) * Param.NxBoxDensity);
  current           = (double***) malloc(sizeof(double**) * Param.NxBoxDensity);
  for(i=0;i<Param.NxBoxDensity;i++){
    //density[i][j] initialized to zero
    density[i]      = (double*) malloc(Param.NyBoxDensity*sizeof(double));
    current[i]      = (double**) malloc(sizeof(double*)*Param.NyBoxDensity);
    for(j=0;j<Param.NyBoxDensity;j++){
      density[i][j]=0;
      current[i][j] = (double*) malloc(2*sizeof(double));
      current[i][j][0]=0;
      current[i][j][1]=0;
    }
  }
  densitycount = 0;
  
  Param.NxBox  = (long) (floor(Param.Lx/Param.rbox)+EPS);
  if(fabs((double)Param.NxBox*Param.rbox-Param.Lx)>EPS){
    printf("Lx has to be a multiple of rbox");
    exit(1);
  }
  
  Param.NyBox = (long) (floor(Param.Ly/Param.rbox)+EPS);
  if(fabs((double)Param.NyBox*Param.rbox-Param.Ly)>EPS){
    printf("Ly has to be a multiple of rbox");
    exit(1);
  }
  
  // Box[i][j]=k if k is the first obstacle in the box.
  // Box[i][j]=-1 if the box contains no obstacles
  Box       = (int**) malloc(Param.NxBox*sizeof(int*));
  // Boxes are initialized as empty
  for (i=0;i<Param.NxBox;i++){
    Box[i] = (int*) calloc(Param.NyBox,sizeof(int));
    for (j=0;j<Param.NyBox;j++)
      Box[i][j]=-1;
  }

  //This is a list of all the boxes around box i in which obstacles
  //may interact with those in box i. epsilonx and epsilony are offset
  //that are used for periodic boundary conditions
  NeighbouringBoxes = (box***) malloc(Param.NxBox*sizeof(box**));
  for (i=0;i<Param.NxBox;i++){
    NeighbouringBoxes[i] = (box**) malloc(Param.NyBox*sizeof(box*));
    for (j=0;j<Param.NyBox;j++){
      NeighbouringBoxes[i][j] = (box*) malloc(9*sizeof(box));
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

      // The box below
      NeighbouringBoxes[i][j][5].i=i;
      // (bool)?(a):(b) is equal to a if bool is true and b otherwise
      NeighbouringBoxes[i][j][5].j= (j==0)?(Param.NyBox-1):(j-1);
      NeighbouringBoxes[i][j][5].epsilonx=0;
      NeighbouringBoxes[i][j][5].epsilony=(j==0)?(-Param.Ly):0;
      
      // The box below, to the left
      NeighbouringBoxes[i][j][6].i= (i==0)?(Param.NxBox-1):(i-1);
      NeighbouringBoxes[i][j][6].j= (j==0)?(Param.NyBox-1):(j-1);
      NeighbouringBoxes[i][j][6].epsilonx=(i==0)?(-Param.Lx):0;
      NeighbouringBoxes[i][j][6].epsilony=(j==0)?(-Param.Ly):0;
      
      // The box to the left
      NeighbouringBoxes[i][j][7].i= (i==0)?(Param.NxBox-1):(i-1);
      NeighbouringBoxes[i][j][7].j=j;
      NeighbouringBoxes[i][j][7].epsilonx=(i==0)?(-Param.Lx):0;
      NeighbouringBoxes[i][j][7].epsilony=0;

      // The box above, to the left
      NeighbouringBoxes[i][j][8].i= (i==0)?(Param.NxBox-1):(i-1);
      NeighbouringBoxes[i][j][8].j= (j==Param.NyBox-1)?(0):(j+1);
      NeighbouringBoxes[i][j][8].epsilonx=(i==0)?(-Param.Lx):0;
      NeighbouringBoxes[i][j][8].epsilony=(j==Param.NyBox-1)?(Param.Ly):0;
    }
  }
  
  /* Place the obstacles in space */
  printf("Place obstacles\n");
  Place_Obstacles(Obstacles,Param,Box);
  
  /* Store obstacle positions */
  for(i=0;i<Param.Nobs;i++){
    
    fprintf(outputobs,"%lg\t%lg\t%lg\n",Obstacles[i].x,Obstacles[i].y,Param.sigma);

    if(Obstacles[i].x<Param.sigma)
      fprintf(outputobs,"%lg\t%lg\t%lg\n",Obstacles[i].x+Param.Lx,Obstacles[i].y,Param.sigma);
    if(Obstacles[i].x>Param.Lx-Param.sigma)
      fprintf(outputobs,"%lg\t%lg\t%lg\n",Obstacles[i].x-Param.Lx,Obstacles[i].y,Param.sigma);

    if(Obstacles[i].y<Param.sigma)
      fprintf(outputobs,"%lg\t%lg\t%lg\n",Obstacles[i].x,Obstacles[i].y+Param.Ly,Param.sigma);
    if(Obstacles[i].y>Param.Ly-Param.sigma)
      fprintf(outputobs,"%lg\t%lg\t%lg\n",Obstacles[i].x,Obstacles[i].y-Param.Ly,Param.sigma);
}
  
  printf("Place particles\n");
  /* Place the run-and-tumble particles in space */
  for(i=0;i<Param.N;i++){
    Place_Particle(i,Particles,Param,Obstacles);
  }
  
  /* Store parameters for reproducibility */
  fprintf(outputparam,"%s\n",command_base);
  for(i=0;i<argc;i++){
    fprintf(outputparam,"%s ",argv[i]);
  }
  fprintf(outputparam,"\n");
  
  fprintf(outputparam,"Number of particles %lg\n", Param.N);
  fprintf(outputparam,"seed is %lld\n", seed);
  fprintf(outputparam,"Total time is %lg\n",totaltime);
  fprintf(outputparam,"Lx is %lg\n", Param.Lx);
  fprintf(outputparam,"Ly is %lg\n", Param.Ly);
  fprintf(outputparam,"Tumbling rate is alpha=%lg\n", Param.alpha);
  fprintf(outputparam,"v0 is %lg\n", Param.v0);
  fprintf(outputparam,"dt is %lg\n", Param.dt);
  fprintf(outputparam,"Number of obstacles is %lg\n", Param.Nobs);
  fprintf(outputparam,"Obstacle radius is sigma=%lg\n", Param.sigma);
  
  fprintf(outputparam,"Time before recording of position is EquilibTimePos=%lg\n", EquilibTimePos);
  fprintf(outputparam,"Interval between recording of positions is StoreInterPos=%lg\n", StoreInterPos);
  
  fprintf(outputparam,"Time before recording of density is EquilibTimeDensity=%lg\n", EquilibTimePos);
  fprintf(outputparam,"Interval between recording of density is InterDensity=%lg\n", StoreInterPos);
  
  fprintf(outputparam,"Bin width for density histogram is %lg\n", Param.dr);
  fflush(outputparam);
  
  /* IMPLEMENTATION OF THE DYNAMICS */
  printf("Implement dynamics\n");
  while(_time < totaltime){
    
    /* MOVE THE PARTICLES */
    Move_Particles_RTP(Particles,Param,&_time,Obstacles,Box,NeighbouringBoxes);

        
    /* RECORD THE POSITIONS */
    //If it is time, store positions
    if(_time>NextStorePos-EPS){
      Record_positions(_time,Particles,Param,outputpos);
      NextStorePos += StoreInterPos;
    }
    
    /* RECORD THE DENSITY */
    if(_time>NextDensity-EPS){
      // Record density and current in the arrays
      Measure_density_current(Particles,density,current,Param);
      // Increase the number of times density & current have been measured
      densitycount += 1.0;
      //Set the time for the next measurement
      NextDensity += InterDensity;
      
      //If it is time, store current and density in the file
      if(_time>NextStoreDensity-EPS){
	//Print data in files and set density[i][j] and current[i][j][0/1] to 0
	Record_density_current(density,current,outputdens,Param,densitycount,_time);
	//Reset densitycount to zero
	densitycount=0;
	//Set the time for the next recording
	NextStoreDensity += StoreInterDensity;
      }
    }
    
    
    /* RECORD THE FLUX */

    /* RECORD THE FORCE */

  }
  
  
  free(Particles);
  free(Obstacles);
  for (i=0;i<Param.NxBox;i++){
    free(Box[i]);
  }
  
  free(Box);
  for (i=0;i<Param.NxBox;i++){
    for (j=0;j<Param.NyBox;j++){
      free(NeighbouringBoxes[i][j]);
    }
    free(NeighbouringBoxes[i]);
  }
  free(NeighbouringBoxes);
  fclose(outputpos);
  fclose(outputobs);
  fclose(outputdens);
    
  for(i=0;i<Param.NxBoxDensity;i++){
    free(density[i]);
    for(j=0;j<Param.NyBoxDensity;j++){
      free(current[i][j]);
    }
    free(current[i]);
  }
  free(density);
  free(current);
  

  return 1;
}















