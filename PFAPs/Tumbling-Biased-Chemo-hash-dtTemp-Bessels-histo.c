/*
  2019-07-11 
  This code simulates the dynamics of run-and-tumble particles whose dynamics is given by
  
  \dot x = v(K \star \rho) \cos\theta
  \dot y = v(K \star \rho) \sin\theta

  at rate alpha, theta becomes theta', where theta' is sampled uniformly in [0,2 pi]

  We use exponential function for v(u) = v_0 \exp(\beta u)

  The kernel K(r)=a(r^2-rmax^2)(r^2-rmin^2). 
  ATTENTION: ICI ON NE RENORMALISE PAS LE KERNEL. 
*/

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
  double time; // time of next tumble
  int bi; // position of box along x axis
  int bj; // position of box along y axis
} particle;

typedef struct box{
  int i; // position of box along x axis
  int j; // position of box along y axis
  double epsilonx; //distance to add along x to take into account PBC
  double epsilony; //distance to add along y to take into account PBC
} box;

// Structure param contains the parameters of the code that will be passed to functions in a compact way
struct param {
  double alpha0; // constant part of the tumbling rate
  double alpha1; //coefficient that modulates the amplitude of the tactic effect on tumbling
  long N;  // number of particles
  double Lx; // system size
  double Ly;
  double dt; // time-step
  double v0; // particle speed
 

  // Parameters of the coarsening Kernel
  double s1; // sqrt(lambda1/nu1) inverse range of first kernel
  double s2; // sqrt(lambda2/nu2) inverse range of second kernel
  double a1; // a /[2 pi nu] amplitude of equivalent potential 
  double a2; // a /[2 pi nu] amplitude of equivalent potential
  double rmin; // small decaying length equal to 1/s1
  
  double rmax; // second change of sign (also the range of the interaction)
  double rmax2; // 
  
  // Parameters of the spatial hashing
  double rbox; //Box size used for hashing. Equal rmax, a priori.
  int NxBox; // number of boxes along x
  int NyBox; // number of boxes along y
};

#define PBC
#define EPS 1e-10

#include "Tumbling-Biased-Chemo-hash-dtTemp-Bessels-histo-functions.c"

int main(int argc, char* argv[]){
  
  /* Variable declaration */
  int i,j,k; // counter
  time_t time_clock; // current physical time

  double _time; // current time
  double rho0;
  double nlocal; // number of particle in a box
  double TimeHisto; // Time after which the histogram is made
  double HistoInter; // Interval between two storage of data for the histogram
  double NextStoreHisto; // Time at which the next data for the histogram are stored.
  double FinalTime; //final time of the simulation
  double StoreInter; // interval at which data are stored
  double NextStore; // next time at which data
  particle* Particles; // arrays containing the particles
  struct param Param; // Structure containing the parameters
  long long seed; // see of the random number generator
  FILE* output, *outputparam,*outputhisto; // Files in which data and parameters are stored
  char name[200]; // string in which the file names are written
  double* GradKer;        // array which contains the gradient of the chemotaxis agent concentration                        // GradKer[2*i] and GradKer[2*i+1] are x and y components acting on particle i
  double* density;
  double* alpha;

  int** Box;            // array containing the first particles in the box
                        // Box[i][j]=k means the first particle in box i,j is k.
  int* Neighbours;   // array containing the neighbours of each
		     // particle. Neighbours[2*i+1]=k means the next
		     // particle after i is k. k=-1 means i is the
		     // last particle. Neighbours[2*i]=k means the
		     // particle before i is k. k=-1 means the particle
		     // is the first in the box.
  box *** NeighbouringBoxes;
  double* dt_temp;// intervalle infinitesimal temporel dependant de la vitesse
  double* NegTumbRateCount; //counter for the percentage of particles having a negative tumbling rate.
  
  /* Read parameters of the simulation*/
  if(argc!=17){
    printf("usage: %s file FinalTime alpha0 alpha1 rho0 Lx Ly dt v0 seed StoreInter s1 s2 a1 a2 TimeHisto\n",argv[0]);//alp[ha and v0 don't appear in the passive code
    exit(1);
  }
  
  i=1;
  
  sprintf(name,"%s-param",argv[i]);
  outputparam=fopen(name,"w");
  sprintf(name,"%s-data",argv[i]);
  output=fopen(name,"w");
  sprintf(name,"%s-histo",argv[i]);
  outputhisto=fopen(name,"w");
  
  i++;
  FinalTime      = strtod(argv[i], NULL); i++;
  Param.alpha0   = strtod(argv[i], NULL); i++;
  Param.alpha1   = strtod(argv[i], NULL); i++;
  rho0           = strtod(argv[i], NULL); i++;  
  Param.Lx       = strtod(argv[i], NULL); i++;
  Param.Ly       = strtod(argv[i], NULL); i++;
  Param.dt       = strtod(argv[i], NULL); i++;
  Param.v0       = strtod(argv[i], NULL); i++;
  seed           = (long long) strtod(argv[i], NULL); i++;
  StoreInter     = strtod(argv[i], NULL); i++;
  Param.s1       = strtod(argv[i], NULL); i++;
  Param.s2       = strtod(argv[i], NULL); i++;
  Param.a1       = strtod(argv[i], NULL); i++;
  Param.a2       = strtod(argv[i], NULL); i++;
  TimeHisto      = strtod(argv[i], NULL); i++;
  
  /* Initialize variables */
  init_genrand64(seed);
  _time               = 0;
  NextStore           = StoreInter;
  HistoInter          = 1;
  NextStoreHisto      = TimeHisto;
  Param.rmax          = 1.;
  Param.rmax2         = Param.rmax*Param.rmax;
  Param.N             = (long) (rho0*Param.Lx*Param.Ly);
  Param.rmin          = 1./Param.s1;
  
  if(1./Param.s2 < 1./Param.s1){
    printf("Use s1>s2\n");
    exit(1);
  }
  
  Particles           = (particle*) malloc(Param.N*sizeof(particle));
  GradKer             = (double*) calloc(2*Param.N,sizeof(double));
  density             = (double*) malloc(Param.N*sizeof(double));
  alpha               = (double*) malloc(Param.N*sizeof(double));
  dt_temp             = (double*) malloc(sizeof(double));
  dt_temp[0]          = Param.dt;
  NegTumbRateCount    = (double*) malloc(sizeof(double));
  NegTumbRateCount[0] = 0;

  // Create box in which to store particle for the spatial hashing
  Param.rbox=Param.rmax;
  Param.NxBox = (int) (floor(Param.Lx/Param.rbox)+EPS);
  if(fabs((double)Param.NxBox*Param.rbox-Param.Lx)>EPS){
    printf("Lx has to be a multiple of rmax");
    exit(1);
  }
  
  Param.NyBox = (int) (floor(Param.Ly/Param.rbox)+EPS);
  if(fabs((double)Param.NyBox*Param.rbox-Param.Ly)>EPS){
    printf("Ly has to be a multiple of rmax");
    exit(1);
  }
  
  Box       = (int**) malloc(Param.NxBox*sizeof(int*));
  for (i=0;i<Param.NxBox;i++){
    Box[i] = (int*) calloc(Param.NyBox,sizeof(int));
    for (j=0;j<Param.NyBox;j++)
      Box[i][j]=-1;
  }
  
  Neighbours = (int*) calloc((Param.N+1)*2,sizeof(int));
  for(i=0;i<Param.N;i++){
    Neighbours[2*i]=-1;
    Neighbours[2*i+1]=-1;
  }
  
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
  for(i=0;i<Param.N;i++){
    //int d=(int) sqrt(Param.N);
    double mindist2,distance2;

    /*I.C. uniformly and randomly drawn*/
    //Particles[i].x     = Param.Lx*genrand64_real3();
    //Particles[i].y     = Param.Ly*genrand64_real3();
    
    /*I.C. : here the particles are disposed according to a grid. But It's really regular only if N is a square !*/
    //Particles[i].x=(i/d)*Param.Lx/d+Param.rmin*genrand64_real3();
    //Particles[i].y=(i%d)*Param.Ly/d+Param.rmin*genrand64_real3();    

    /* Here we uniformly randomly drawn particles one by one, rejecting a particle's initial position if it's too close to another.To be used for hard core interaction.*/
    mindist2 = 0;//initialization so that we pass at least once in the loop.
    while(mindist2<Param.rmin/20)//Critère un peu arbitraire...
      {
	mindist2=Param.Lx*Param.Ly;
	Particles[i].x  = Param.Lx*genrand64_real3();
	Particles[i].y  = Param.Ly*genrand64_real3();
	for(j=0;j<i;j++)
	  {
	    distance2 = (Particles[i].x-Particles[j].x)* (Particles[i].x-Particles[j].x)+(Particles[i].y-Particles[j].y)* (Particles[i].y-Particles[j].y);
	    if(distance2<mindist2)
	      {
		mindist2=distance2;
	      }
	  }
      }   
    
    Particles[i].theta = 2*M_PI*genrand64_real2();

    //Add particle in the good box.
    Particles[i].bi    = floor(Particles[i].x/Param.rbox);
    Particles[i].bj    = floor(Particles[i].y/Param.rbox);
    AddinBox(i,Particles[i].bi,Particles[i].bj,Box,Neighbours);
  }


  
  /* Store parameters */
  fprintf(outputparam,"usage: %s file FinalTime alpha0 alpha1 rho0 Lx Ly dt v0 seed StoreInter s1 s2 a1 a2 TimeHisto\n",argv[0]);
  
  for(i=0;i<argc;i++){
    fprintf(outputparam,"%s ",argv[i]);
  }
  fprintf(outputparam,"\n");
  
  fprintf(outputparam,"FinalTime is %lg\n",FinalTime);
  fprintf(outputparam,"alpha0 is %lg\n",Param.alpha0);
  fprintf(outputparam,"alpha1 is %lg\n",Param.alpha1);
  fprintf(outputparam,"rho0 is %lg\n", rho0);
  fprintf(outputparam,"N is %ld\n", Param.N);
  fprintf(outputparam,"Lx is %lg\n", Param.Lx);
  fprintf(outputparam,"Ly is %lg\n", Param.Ly);
  fprintf(outputparam,"dt is %lg\n", Param.dt);
  fprintf(outputparam,"v0 is %lg\n", Param.v0);
  fprintf(outputparam,"seed is %lld\n", seed);
  fprintf(outputparam,"StoreInter is %lg\n", StoreInter);
  fprintf(outputparam,"rmin is %lg\n", Param.rmin);
  fprintf(outputparam,"rmax is %lg\n", Param.rmax);
  fprintf(outputparam,"s1 is %lg\n", Param.s1);
  fprintf(outputparam,"s2 is %lg\n", Param.s2);
  fprintf(outputparam,"a1 is %lg\n", Param.a1);
  fprintf(outputparam,"a2 is %lg\n", Param.a2);
  
  fprintf(outputparam,"TimeHisto is %lg\n", TimeHisto);  
  fprintf(outputparam,"packing fraction large size is %lg\n", Param.N*M_PI/Param.s2/Param.s2/Param.Lx/Param.Ly/4);
  fprintf(outputparam,"packing fraction small size is %lg\n", Param.N*M_PI/Param.s1/Param.s1/Param.Lx/Param.Ly/4);
  
  fflush(outputparam);

  // Save initial condition
  /*  Density_TopHat_Hashing_PBC(Particles,Param,density,Box,Neighbours,NeighbouringBoxes); */
  /*  Kernel_TwoSchwartz_Hashing_PBC(Particles,Param,velocity,Box,Neighbours,NeighbouringBoxes); */
  
  /* for(i=0;i<Param.N;i++){ */
  /*   fprintf(output,"%lg\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\n",_time,i,Particles[i].x,Particles[i].y,Particles[i].theta,density[i],Param.v0*exp(Param.beta*velocity[i])); */
  /* } */
  /* fprintf(output,"\n"); */
  /* fflush(output); */

  time_clock       = time(NULL);
  
  /* Run the dynamics */
  
  while(_time<FinalTime){
    // Move the particles
    Move_Particles_RTP(Particles,Param,GradKer,alpha,_time,Box,Neighbours,NeighbouringBoxes,dt_temp,NegTumbRateCount);
    _time += dt_temp[0];
    
    //If it is time, store data
    if(_time>NextStore-EPS){
      Density_TopHat_Hashing_PBC(Particles,Param,density,Box,Neighbours,NeighbouringBoxes);
      // Kernel_TwoSchwartzShoulders_Hashing_PBC(Particles,Param,velocity,Box,Neighbours,NeighbouringBoxes);
	
      // Store positions, angles, density and speeds
      for(i=0;i<Param.N;i++){
	fprintf(output,"%lg\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\n",_time,i,Particles[i].x,Particles[i].y,Particles[i].theta,density[i],alpha[i]);
      }
      fprintf(output,"\n");
      fflush(output);
      //Store dt_temp[0] and NegTumbRateCount[0]
      fprintf(outputparam,"#time: %lg dt: %lg\t Percentage of particles having negative tumbling rate at this time: %lg\n",_time,dt_temp[0],100*NegTumbRateCount[0]/Param.N);
      fflush(outputparam);
      
      //Increase NextStore
      NextStore += StoreInter;
    }

    
    //If it's time, record the number of particles in each box, rescaled by the box'size
    if(_time>NextStoreHisto-EPS){
      for (i=0;i<Param.NxBox;i++)
	for(j=0;j<Param.NyBox;j++){
	  nlocal=0; // local number of particles
	  k=Box[i][j]; // first particle in the box
	  // as long as this is a particle, increase nlocal and move down the list
	  while(k!=-1){
	    nlocal++;
	    k=Neighbours[2*k+1];
	  }
	  fprintf(outputhisto,"%lg\t%lg\n",_time,nlocal/Param.rbox/Param.rbox);
	}
      NextStoreHisto+=HistoInter;
      fprintf(outputhisto,"#\n");
    }
    
  }
  
  printf("Simulation time: %ld seconds\n",time(NULL)-time_clock);
  fprintf(outputparam,"#Simulation time: %ld seconds\n",time(NULL)-time_clock);
  
  free(Neighbours); 
  for (i=0;i<Param.NxBox;i++){
    free(Box[i]);
  }
  free(Box);
  for (i=0;i<Param.NxBox;i++){
    for (j=0;j<Param.NyBox;j++)
      free(NeighbouringBoxes[i][j]);
    free(NeighbouringBoxes[i]);
  }
  free(NeighbouringBoxes);
  free(Particles);
  free(GradKer);
  free(density);
  free(dt_temp);
  free(NegTumbRateCount);
 

  return 1;
}


