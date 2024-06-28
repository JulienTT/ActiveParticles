void ERROR(char* msg){
  printf(msg);
  exit(1);
}

#if defined(WCA) || defined(LJ)
#include "WCA.c"
#endif

/*
  To add someone in a box, you need to put it there and to set its new
  neighbours. This functions add at the top.
*/
void AddinBox(long i,long bi, long bj,long** Box,long* Neighbours){
  long k;
  // Save the old state of Box, -1 or a particle number
  k=Box[bi][bj];
  
  //Put particle i at the top
  Box[bi][bj]=i;
  // There is no one above particle i
  Neighbours[2*i]=-1;
  // k is after i
  Neighbours[2*i+1]=k;
  // If k is particle, i is before k
  if(k!=-1)
    Neighbours[2*k]=i;
}

/* 
   To remove someone from the box, you have to remove it, but you do
   not need to update its neighbours. This will be done when it is
   added to a new box
*/
void RemovefromBox(long i,long bi,long bj,long** Box,long* Neighbours){
  long next;
  long prev;
  // Store the next particle
  next=Neighbours[2*i+1];
  // If the particle is the first in the box
  if(Box[bi][bj]==i){
    // next is the new top of the box
    Box[bi][bj]=next;
    // if next is a particle, there is no one before it
    if(next!=-1){
      Neighbours[2*next]=-1;
    }
  }
  // If there is someone before you
  else{
    // Store the previous particle
    prev=Neighbours[2*i];
    // The next of the previous is your next
    Neighbours[2*prev+1]=next;
    // If next is a particle, its previous is your previous
    if(next!=-1)
      Neighbours[2*next]=prev;
  }
}

void Density_TopHat_Hashing_PBC(particle* Particles,param Param,double* density,long** Box,long* Neighbours,box *** NeighbouringBoxes){
  long i,j,k,bi,bj,nbi,nbj;
  double dist2;
  double dx,dy,epsilonx,epsilony;
  
  //Reset the array density at [0,0,0,...,0]
  memset(density, 0, sizeof(double)*Param.N);
  
  //Cycle through all boxes
  for(bi=0;bi<Param.NxBox;bi++){
    for(bj=0;bj<Param.NyBox;bj++){
      
      //Take the first particle in the box
      i=Box[bi][bj];
      
      // As long as there are particle in the box, update their
      // densities and their contribution to others
      while(i!=-1){
	
	/* First look at the current box */
	j=Neighbours[2*i+1];
	
	// While particle j has neighbours
	while(j!=-1){
	  dx=fabs(Particles[i].x-Particles[j].x);
	  dy=fabs(Particles[i].y-Particles[j].y);
	  dist2=dx*dx+dy*dy;
	  // if j contributes to density[i], the converse is also true
	  // Here we use a top-hat kernel
	  if(dist2<Param.rbox2){
	    density[i] += 1.;
	    density[j] += 1.;
	  }
	  // Now that j has contributed, look at the next particle
	  j=Neighbours[2*j+1];
	}
	
	/* cycle through the 4 neighbouring boxes to compute their contributions */
	
	for(k=1;k<5;k++){
	  // Store the box 
	  nbi=NeighbouringBoxes[bi][bj][k].i;
	  nbj=NeighbouringBoxes[bi][bj][k].j;
	  epsilonx=NeighbouringBoxes[bi][bj][k].epsilonx;
	  epsilony=NeighbouringBoxes[bi][bj][k].epsilony;
	  
	  // Take the first particle of the neighbouring box
	  j=Box[nbi][nbj];
	  //Cycle through all particles after j
	  while(j!=-1){
	    dx=fabs(Particles[j].x+epsilonx-Particles[i].x);
	    dy=fabs(Particles[j].y+epsilony-Particles[i].y);
	    dist2=dx*dx+dy*dy;
	    // if j contributes to density[i], the converse is also true
	    // Here we use a top-hat kernel
	    if(dist2<Param.rbox2){
	      density[i] += 1.;
	      density[j] += 1.;
	    }
	    j=Neighbours[2*j+1];
	  }
	}
	i=Neighbours[2*i+1];
      }
    }
  }
  for(i=0;i<Param.N;i++)
    density[i] = density[i]/(M_PI*Param.rbox2);
}



void Move_Particles_ABP(particle* Particles,param Param,double* forces,double _time,long** Box,long* Neighbours,box *** NeighbouringBoxes, FILE* outputparam
#ifdef PRESSURE
			,double* pressure,double* average_density	
#endif
#ifdef PRESSURE
			,double* pressure,double* average_density	
#endif
#ifdef STRESSTENSOR
			,double*** sigma_IK,double** sigma_a,double* next_measure_stress_tensor,double* next_store_stress_tensor,double inter_store_stress_tensor,long* count_measure_stress_tensor
#endif
#ifdef OBSERVABLE
			,double* Energy	
#endif
			){
  long i;
  double newbi, newbj; // New boxes where particles will be moved
  double dr2,dr2_max;  // Used to compute the largest displacement
  double dx,dy;        // Used to compute displacements
  long i_max;          // index of particle with largest displacement
  

#ifdef STRESSTENSOR
  int StressBool=0;
  if(_time>next_measure_stress_tensor){
    StressBool=1;
    next_measure_stress_tensor[0]+=inter_store_stress_tensor;
  }
#endif

#if defined(WCA) || defined(LJ)
  /*Compute force for each particle*/
  Force_LJ(Particles,Param,forces,Box,Neighbours,NeighbouringBoxes
#ifdef STRESSTENSOR
	   ,StressBool,sigma_IK,sigma_a,count_measure_stress_tensor
#endif
#ifdef OBSERVABLE
	   ,Energy
#endif
	   );

#endif
#endif
  
  /*Largest displacement*/
  dr2_max = 0;
  
  for(i=0;i<Param.N;i++){
    //Compute displacements along x and y
    
    dx = Param.v0 * cos(Particles[i].theta) * Param.dt + Param.sqrt2Dtdt*gasdev();
    dy = Param.v0 * sin(Particles[i].theta) * Param.dt + Param.sqrt2Dtdt*gasdev();
    
#if defined(LJ)|| defined(WCA)
    //printf("force applied\n");
    dx += Param.mu * forces[2*i] * Param.dt;
    dy += Param.mu * forces[2*i+1] * Param.dt;
#endif
    
#ifdef CLOSEDBC
    if(Particles[i].x>Param.x_wall_right){
      //Absolute value of the force
      double fwalldt= Param.Omega * pow( Particles[i].x-Param.x_wall_right ,Param.nu-1.) * Param.dt;
      dx          -= Param.mu * fwalldt;
      pressure[1] += fwalldt;
    }
    if(Particles[i].x<Param.x_wall_left){
      //Absolute value of the force
      double fwalldt= Param.Omega * pow( Param.x_wall_left-Particles[i].x  ,Param.nu-1.) * Param.dt;
      dx          += Param.mu * fwalldt;
      pressure[0] += fwalldt;
    }
#endif
    
    // squared displacement
    dr2=dx*dx+dy*dy;
    
    // update largest displacement
    if (dr2>dr2_max)
      dr2_max=dr2;
    
    Particles[i].x     += dx;
    Particles[i].y     += dy;
    Particles[i].theta += Param.sqrt2Drdt*gasdev();
    
    if(Particles[i].theta>2*M_PI)
      Particles[i].theta-=2*M_PI;
    if(Particles[i].theta<0)
      Particles[i].theta+=2*M_PI;
    
    //Correct displacement for periodic boundary conditions
#ifdef PBC
    if(Particles[i].x>Param.Lx){
      Particles[i].x -= Param.Lx;
      if(Particles[i].x>Param.Lx){
	fprintf(outputparam,"x>2Lx at time %lg\n",_time);
	ERROR("x>2Lx\n");
      }
    }
    if(Particles[i].x<0){
      Particles[i].x += Param.Lx;
      if(Particles[i].x<0){
	fprintf(outputparam,"x<-Lx at time %lg\n",_time);
	ERROR("x<-Lx\n");
      }
    }
#endif
    //If closed, check you are not too close to the PBC
#ifdef CLOSEDBC
    if(Particles[i].x>Param.Lx-Param.rbox){
      fprintf(outputparam,"x_wall_right too small\n");
      ERROR("x_wall_right too large\n");
    }
    if(Particles[i].x<Param.rbox){
      fprintf(outputparam,"x_wall_left too small\n");
      ERROR("x_wall_left too small\n");
    }
#endif
    //Always use PBC along y
    if(Particles[i].y>Param.Ly){
      Particles[i].y -= Param.Ly;
      if(Particles[i].y>Param.Ly){
	fprintf(outputparam,"y>2Ly at time %lg\n",_time);
	ERROR("y>2Ly\n");
      }
    }
    if(Particles[i].y<0){
      Particles[i].y += Param.Ly;
      if(Particles[i].y<0){
	fprintf(outputparam,"y<-Ly at time %lg\n",_time);
	ERROR("y<-Ly\n");
      }
    }
    
#ifdef PRESSURE
    if( Particles[i].x>Param.Lx*2./5. && Particles[i].x<=Param.Lx*3./5. )
      average_density[0] += Param.dt;
#endif
    
    //Test if the particle has changed box
    newbi    = floor(Particles[i].x/Param.rbox);
    newbj    = floor(Particles[i].y/Param.rbox);
    
    if(Particles[i].bi!=newbi||Particles[i].bj!=newbj){
      RemovefromBox(i,Particles[i].bi,Particles[i].bj,Box,Neighbours);
      AddinBox(i,newbi,newbj,Box,Neighbours);
      Particles[i].bi=newbi;
      Particles[i].bj=newbj;
    }
  }
  
#if defined(WCA) || defined(LJ)
  //If the largest displacement is larger than sigma/10, say so
  if(dr2_max>Param.sigma2*0.01)
    fprintf(outputparam,"at time %lg, dr_max=%lg largest acceptable:%lg\n",_time,sqrt(dr2_max),Param.sigma*0.1);
#endif
}

double Distance2(particle* Particles,long i,long j, param Param){
  double dx=fabs(Particles[i].x-Particles[j].x);
  if(dx>Param.Lx/2.)
    dx=Param.Lx-dx;
  double dy=fabs(Particles[i].y-Particles[j].y);
  if(dy>Param.Ly/2.)
    dy=Param.Ly-dy;
  return dx*dx+dy*dy;
}
 
double SmallestDistance(particle* Particles,param Param){
  long i,j;
  double rmin2=Param.Lx*Param.Lx;
  double dx,dy,dr2;
  
  for(i=0;i<Param.N;i++)
    for(j=i+1;j<Param.N;j++){
      dr2=Distance2(Particles,i,j,Param);
      if(dr2<rmin2)
	rmin2=dr2;
    }
  return sqrt(rmin2);

}

/*Check input parameters, attribute them and store them*/
void  CheckAndReadInput(char* command_base,int argc, char* argv[], double* rho0, param* Param,
#ifdef _MT
			long long* seed,
#endif
#ifdef _PCG
			pcg128_t* seed,
#endif
			double* FinalTime, double* EquilibTimePos, double* StoreInterPos, int* width, double* rhomax, double* HistoInter, double* StoreHistoInter, FILE** outputparam, FILE** outputpos, FILE** outputhisto
#ifdef PRESSURE
			, FILE** outputpressure, double* StorePressureInter, double* average_density
#endif
#ifdef STRESSTENSOR
			, FILE** output_stress_tensor, double* Inter_Measure_Stress_Tensor,double* Inter_Store_Stress_Tensor
#endif
#ifdef GIVENIC
			, FILE** input
#endif
#ifdef SLABIC
			, double* rhog, double* rhol, double* liquidfraction
#endif
#ifdef OBSERVABLE
			, FILE** outputenergy, double* StoreInterEnergy
#endif
			
			){
  int i;
  int argctarget=21;
  char name[200];      // string in which the file names are written
  strcat(command_base, "usage: ");
  strcat(command_base, argv[0]);
  strcat(command_base," file rho0 Lx Ly dt v0 Dr Dt mu sigma epsilon seed FinalTime EquilibTimePos StoreInterPos rbox width rhomax HistoInter StoreHistoInter");
  
#ifdef GIVENIC
  strcat(command_base," input N");
  argctarget += 2;
#endif
  
#ifdef SLABIC
  strcat(command_base," rho_g rho_ell liquid_fraction");
  argctarget += 3;
#endif
  
#ifdef CLOSEDBC
  strcat(command_base," Omega nu x_wall_left x_wall_right");
  argctarget += 4;
#endif
  
#ifdef PRESSURE
  strcat(command_base," StorePressureInter");
  argctarget += 1;
#endif

#ifdef STRESSTENSOR
  strcat(command_base," Inter_Measure_Stress_Tensor Inter_Store_Stress_Tensor");
  argctarget += 2;
#endif

#ifdef OBSERVABLE
  strcat(command_base," StoreInterEnergy");
  argctarget += 1;
#endif
  
  if(argc!=argctarget){
    printf("%s\n",command_base);
    exit(1);
  }
  
  i=1;
  //File where parameters are stored
  sprintf(name,"%s-param",argv[i]);
  outputparam[0]=fopen(name,"w");
  
  //File where positions are stored
  sprintf(name,"%s-pos",argv[i]);
  outputpos[0]=fopen(name,"w");

  //File where histogram is stored
  sprintf(name,"%s-histo",argv[i]);
  outputhisto[0]=fopen(name,"w");
  
#ifdef PRESSURE
  //File where pressure is stored
  sprintf(name,"%s-pressure",argv[i]);
  outputpressure[0]=fopen(name,"w");
#endif

#ifdef STRESSTENSOR
  //File where pressure is stored
  sprintf(name,"%s-stress-tensor",argv[i]);
  output_stress_tensor[0]=fopen(name,"w");
#endif

#ifdef OBSERVABLE
  //File where potential energy is stored
  sprintf(name,"%s-energy",argv[i]);
  outputenergy[0]=fopen(name,"w");
#endif
  
  i++;
  rho0[0]              = strtod(argv[i], NULL); i++;  
  Param[0].Lx          = strtod(argv[i], NULL); i++;
  Param[0].Ly          = strtod(argv[i], NULL); i++;
  Param[0].dt          = strtod(argv[i], NULL); i++;
  Param[0].v0          = strtod(argv[i], NULL); i++;
  Param[0].Dr          = strtod(argv[i], NULL); i++;
  Param[0].Dt          = strtod(argv[i], NULL); i++;
  Param[0].mu          = strtod(argv[i], NULL); i++;
  Param[0].sigma       = strtod(argv[i], NULL); i++;
  Param[0].epsilon     = strtod(argv[i], NULL); i++; 
#ifdef _MT
  seed[0]              = (long long) strtod(argv[i], NULL); i++;
#endif
#ifdef _PCG
  seed[0]              = (pcg128_t) strtod(argv[i], NULL); i++;
#endif
  
  FinalTime[0]         = strtod(argv[i], NULL); i++;
  EquilibTimePos[0]    = strtod(argv[i], NULL); i++;
  StoreInterPos[0]     = strtod(argv[i], NULL); i++;
  Param[0].rbox        = strtod(argv[i], NULL); i++;
  width[0]             = (int) strtod(argv[i], NULL); i++;
  rhomax[0]            = strtod(argv[i], NULL); i++;
  HistoInter[0]        = strtod(argv[i], NULL); i++;
  StoreHistoInter[0]   = strtod(argv[i], NULL); i++;
  
#ifdef GIVENIC
  sprintf(name,"%s",argv[i]);i++;
  printf("read from file %s\n",name);
  Param[0].N           = (long) strtod(argv[i], NULL); i++;
  printf("N is %ld\t rho is %lg\n",Param[0].N,Param[0].N/Param[0].Lx/Param[0].Ly);
  input                = fopen(name,"r");
#endif
  
#ifdef SLABIC
  rhog[0]              = strtod(argv[i], NULL); i++;
  rhol[0]              = strtod(argv[i], NULL); i++;
  liquidfraction[0]    = strtod(argv[i], NULL); i++;
#endif

#ifdef CLOSEDBC
  Param[0].Omega       = strtod(argv[i], NULL); i++;
  Param[0].nu          = strtod(argv[i], NULL); i++;
  Param[0].x_wall_left = strtod(argv[i], NULL); i++;
  Param[0].x_wall_right= strtod(argv[i], NULL); i++;
#endif
  
#ifdef PRESSURE
  StorePressureInter[0]= strtod(argv[i], NULL); i++;
  average_density[0]   = 0;
#endif

#ifdef STRESSTENSOR
  Inter_Measure_Stress_Tensor[0] = strtod(argv[i], NULL); i++;
  Inter_Store_Stress_Tensor[0] = strtod(argv[i], NULL); i++;
#endif
  
#ifdef OBSERVABLE
  StoreInterEnergy[0]  = strtod(argv[i], NULL); i++;
#endif
  
}


void DefineNeighbouringBoxes(box *** NeighbouringBoxes,param Param){
  int i,j;
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
}


/* 
   This function randomly draws particles one by one, rejecting a
   particle's initial position if it's too close to another. Can be used
   for hard core interaction at low and moderate densities.
*/
void Place_Random_IC_Continuous(param Param,particle* Particles,long** Box,long* Neighbours){
  long i,j;
  double mindist2,distance2;

  for(i=0;i<Param.N;i++){
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
    //printf("x\t%lg\ty\t%lg\n",Particles[i].x,Particles[i].y);
    Particles[i].theta = 2*M_PI*genrand64_real2();
    
    //Add particle in the good box.
    Particles[i].bi    = floor(Particles[i].x/Param.rbox);
    Particles[i].bj    = floor(Particles[i].y/Param.rbox);
    AddinBox(i,Particles[i].bi,Particles[i].bj,Box,Neighbours);
  }
  printf("Smallest distance %lg\n",SmallestDistance(Particles,Param));
}

/*
  This function reads the initial condition in the file input.
*/
void Place_Given_IC(param Param,FILE* input,particle* Particles,long** Box,long* Neighbours){
  long i;
  double x,y,theta;
  for(i=0;i<Param.N;i++){
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
    else
      ERROR("Not the right number of lines in input file\n");
  }
}

/*
  We create two hexagonal lattices of widths (Lx liquidfraction) and
  Lx (1-liquidfraction) to place liquid and gas particles.
  
  The lattice spacing is 0.9 sigma so that particles do not interact
  too strongly.
  
  The number of particles in the phases are Nliquid and Ngas,
  respectively
*/
void Place_Slab_IC_Lattice(param Param,double liquidfraction,FILE* outputparam,particle* Particles,long** Box,long* Neighbours,long Ngas, long Nliquid){
  double a; //horizontal lattice spacing
  double h; //heiht between two layers. Equal to a*cos(pi/6)
  int NxL; //number of rows in the lattice for the liquid phase
  int NxG; //number of rows in the lattice for the gas phase
  int Ny;  //number of lines in the lattices
  long NsitesGas; //Total number of available sites in the gas phase
  long NsitesLiq; //Total number of available sites in the liquid phase
  int* GasPhase; //Lattice of gas phase
  int* LiquidPhase;//Lattice of gas phase
  long nb; //Number of particles placed so far
  long i,j,k;
  
  a = Param.sigma;
  h = a*cos(M_PI/6.);
  //Number of sites in x direction in liquid phase
  NxL = (floor) (Param.Lx * liquidfraction     / a );
  //Number of sites in x direction in gas phase
  NxG = (floor) (Param.Lx * (1-liquidfraction) / a );
  //Number of sites in y direction
  Ny  = (floor) (Param.Ly / h );
  //Total number of free sites in each phase
  NsitesGas=NxG*Ny;
  NsitesLiq=NxL*Ny;
  //We create the corresponding lattices, which are currently empty
  GasPhase    = (int*) calloc(NsitesGas,sizeof(int));
  LiquidPhase = (int*) calloc(NsitesLiq,sizeof(int));
  //if there are enough spaces for the Ngas particles
  if(Ngas<NsitesGas){
    for(i=0;i<Ngas;i++){
      int bool=1;
      // We pull a site at random and place the particle there if the
      // site is not occupied
      while(bool==1){
	k=genrand64_int64()%NsitesGas;
	if(GasPhase[k]==0){
	  GasPhase[k]=1;
	  bool=0;
	}
      }
    }
  }
  else
    ERROR("Not enough gas sites\n");

  //if there are enough sites for the Nliquid particles
  if(Nliquid<NsitesLiq){
    for(i=0;i<Nliquid;i++){
      int bool=1;
      // We pull a site at random and place the particle there if the
      // site is not occupied
      while(bool==1){
	k=genrand64_int64()%NsitesLiq;
	if(LiquidPhase[k]==0){
	  LiquidPhase[k]=1;
	  bool=0;
	}
      }
    }
  }
  else
    ERROR("Not enough liquid sites\n");
  
  fprintf(outputparam,"Actual gas density is %lg\n", Ngas/Param.Lx/Param.Ly*(1-liquidfraction));
  fprintf(outputparam,"Actual liquid density  is %lg\n", Nliquid/Param.Lx/Param.Ly*(1-liquidfraction));
  
  // We now place particles in Particles, given the arrays
  //nb is the label of the particle to be placed. 
  nb=0;
  //We start with the particles in NsitesLiq. We go through all the
  //sites of the lattice LiquidPhase.
  for(i=0;i<NsitesLiq;i++){
    //If it is occupied, place the particle.
    if(LiquidPhase[i]==1){
      //we know that i = k + NxL * j
      k=i%NxL;
      j=(i-k)/NxL;
      //Compute the corresponding x. The left end of the liquid phase
      //is at Lx (1-liquidfraction). Depending on whether the line is
      //even or odd, there is an a/2 offset
      Particles[nb].x= (1-liquidfraction) * Param.Lx * .5 + k*a + a * .5 * (j%2) ;
      Particles[nb].y=j*h;
      Particles[nb].theta = 2*M_PI*genrand64_real2();

      //Add particle in the good box.
      Particles[nb].bi    = floor(Particles[nb].x/Param.rbox);
      Particles[nb].bj    = floor(Particles[nb].y/Param.rbox);
      AddinBox(nb,Particles[nb].bi,Particles[nb].bj,Box,Neighbours);
      nb++;
    }
  }
  //Let's do the same with the gas phase. x starts at (1+liquid
  //fraction)Lx but we have to use periodic boundary conditions
  for(i=0;i<NsitesGas;i++){
    if(GasPhase[i]==1){
      k=i%NxG;
      j=(i-k)/NxG;
      Particles[nb].x= (1+liquidfraction) * Param.Lx *.5 + k*a + a * .5 * (j%2);
      if(Particles[nb].x>Param.Lx)
	Particles[nb].x -= Param.Lx;
      Particles[nb].y=j*h;
      Particles[nb].theta = 2*M_PI*genrand64_real2();

      //Add particle in the good box.
      Particles[nb].bi    = floor(Particles[nb].x/Param.rbox);
      Particles[nb].bj    = floor(Particles[nb].y/Param.rbox);
      AddinBox(nb,Particles[nb].bi,Particles[nb].bj,Box,Neighbours);

      nb++;
    }
  }
  fprintf(outputparam,"Total number of particles placed %d\n",nb);
  fprintf(outputparam,"Planned number of particles %ld\n",Ngas+Nliquid);
  if(nb!=Param.N)
    ERROR("Wrong number of particles placed\n");
  free(GasPhase);
  free(LiquidPhase);
}


void StoreParam(param Param,FILE* outputparam,char* command_base,int argc, char* argv[],double rho0,
#ifdef _MT
		long long seed,
#endif
#ifdef _PCG
		pcg128_t seed,
#endif
		double FinalTime,double EquilibTimePos,double StoreInterPos,int width, double rhomax, double HistoInter,int Nbin,double StoreHistoInter
#ifdef PRESSURE
		, double StorePressureInter
#endif
#ifdef STRESSTENSOR
	     , double Inter_Measure_Stress_Tensor, double Inter_Store_Stress_Tensor
#endif
#ifdef OBSERVABLE
		, double StoreInterEnergy
#endif
		){
  int i;
  
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
  fprintf(outputparam,"Mobility is %lg\n", Param.mu);
  fprintf(outputparam,"sigma is %lg\n", Param.sigma);
  fprintf(outputparam,"epsilon is %lg\n", Param.epsilon);
#ifdef _MT
  fprintf(outputparam,"seed is %lld\n", seed);
#endif
#ifdef _PCG
  fprintf(outputparam,"seed is %" PRIu64 "\n", (uint64_t) seed);
#endif
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

#ifdef STRESSTENSOR
  fprintf(outputparam,"Inter_Measure_Stress_Tensor is %lg\n",Inter_Measure_Stress_Tensor);
  fprintf(outputparam,"Inter_Store_Stress_Tensor is %lg\n",Inter_Store_Stress_Tensor);
#endif

#ifdef OBSERVABLE
  fprintf(outputparam,"StoreInterEnergy is %lg\n",StoreInterEnergy);
#endif
  
  fprintf(outputparam,"nbins is %d\n", Nbin);
#if defined(WCA) || defined(LJ)
  fprintf(outputparam,"rmax2 is %lg\n", Param.rmax2);
#endif
  fprintf(outputparam,"rbox2 is %lg\n", Param.rbox2);
  fprintf(outputparam,"N is %ld\n", Param.N);
  fprintf(outputparam,"normalized packing fraction is N sigma^2 /(Lx Ly) = %lg\n", Param.N*Param.sigma*Param.sigma/Param.Lx/Param.Ly);
  fflush(outputparam);
}

void StorePositions(particle* Particles, param Param, double* density,long ** Box,long* Neighbours, box*** NeighbouringBoxes,FILE* outputpos,double _time,double* NextStorePos,double StoreInterPos){
  long i;
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
  NextStorePos[0] += StoreInterPos;
}

void UpdateHistogram(param Param, int width,long* Neighbours,int Nbin,double* histogram,double* histocount,double* NextHisto,double HistoInter,long** Box,FILE* outputparam,double rhomax,double _time){
  int i,j,i1,j1;
  long n0,part;
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
	histocount[0] += 1;
      }
      else{
	printf("rhomax too small, n0 is %d and Nbin is %d\n",n0,Nbin);
	fprintf(outputparam,"rhomax= %lg too small at time %lg\n",rhomax,_time);
	fflush(outputparam);
	exit(1);
      }
    }
  }
  NextHisto[0] += HistoInter;
}

void StoreHistogram(int Nbin,FILE* outputhisto,double _time,double boxarea,double* histogram,double* histocount, double* NextStoreHisto, double StoreHistoInter,param Param){
  int i;
  for(i=0;i<Nbin;i++){
    fprintf(outputhisto,"%lg\t%lg\t%lg\n",_time,(double) i/boxarea*Param.sigma2,histogram[i]/histocount[0]*boxarea/Param.sigma2);
  }
  NextStoreHisto[0] += StoreHistoInter;
  memset(histogram, 0, sizeof(double)*Nbin);
  histocount[0]=0;
  fprintf(outputhisto,"\n");
  fflush(outputhisto);
}

void StorePressure(FILE* outputpressure, double _time, double* pressure, param Param,double* average_density,double StorePressureInter,double* NextStorePressure){

  fprintf(outputpressure,"%lg\t%lg\t%lg\t%lg\n",_time,pressure[0]/StorePressureInter/Param.Ly,pressure[1]/StorePressureInter/Param.Ly,average_density[0]/Param.Lx*5/Param.Ly/StorePressureInter);
      
  NextStorePressure[0] += StorePressureInter;
  average_density[0]    = 0;
  pressure[0]           = 0;
  pressure[1]           = 0;
  fflush(outputpressure);
}

void StoreEnergy(FILE* outputenergy, double _time, double* Energy, param Param,double StoreInterEnergy,double* NextStoreEnergy){

  fprintf(outputenergy,"%lg\t%lg\n",_time,Energy[0]/StoreInterEnergy/Param.N);
  NextStoreEnergy[0] += StoreInterEnergy;
  Energy[0]           = 0;
  fflush(outputenergy);
}
