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

void Density_TopHat_Hashing_PBC(particle* Particles,struct param Param,double* density,long** Box,long* Neighbours,box *** NeighbouringBoxes){
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

/*
  This function finds dt such that

  orderdt^2 dt^2 + orderrootdt^2 dt = rmax
  
 */
double FindDt(particle* Particles,double* orderdt,double* orderrootdt,long i_max,struct param Param){
  double a,b;
  a = orderdt[2*i_max]*orderdt[2*i_max] + orderdt[2*i_max+1]*orderdt[2*i_max+1];
  b = orderrootdt[2*i_max]*orderrootdt[2*i_max] + orderrootdt[2*i_max+1]*orderrootdt[2*i_max+1];
  return (-b + sqrt(b*b+4.*a*Param.rmax2) ) / (2.*a);
};

void WCA_Force(double* Force,double dx,double dy,double dist2,struct param Param){
  double dist6=dist2*dist2*dist2;
  double amp=Param.amp/(dist6*dist2)*(1-2*Param.sigma6/dist6);
  Force[0]=dx*amp;
  Force[1]=dy*amp;
}

void Force_LJ(particle* Particles,struct param Param, double* forces,long** Box,long* Neighbours,box*** NeighbouringBoxes){
  
  long i,j,k,bi,bj,nbi,nbj;
  double dist2;
  double* Force = (double*) calloc(2,sizeof(double));
  double dx,dy,epsilonx,epsilony;
  
  //Reset the array density at [0,0,0,...,0]
  memset(forces, 0, sizeof(double)*2*Param.N);
  
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
	  dx=Particles[j].x-Particles[i].x;
	  dy=Particles[j].y-Particles[i].y;
	  dist2=dx*dx+dy*dy;
	  // if j contributes to density[i], the converse is also true
	  // Here we use a top-hat kernel
	  if(dist2<Param.rmax2){
	    //This compute the force on particle i
	    WCA_Force(Force,dx,dy,dist2,Param);
	    forces[2*i]   += Force[0];
	    forces[2*i+1] += Force[1];

	    //The force on j is minus the force on i
	    forces[2*j]   -= Force[0];
	    forces[2*j+1] -= Force[1];
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
	    dx=Particles[j].x+epsilonx-Particles[i].x;
	    dy=Particles[j].y+epsilony-Particles[i].y;
	    dist2=dx*dx+dy*dy;
	    // if j contributes to density[i], the converse is also true
	    // Here we use a top-hat kernel
	    if(dist2<Param.rmax2){
	      //This compute the force on particle i
	      WCA_Force(Force,dx,dy,dist2,Param);
	      forces[2*i]   += Force[0];
	      forces[2*i+1] += Force[1];
	      forces[2*j]   -= Force[0];
	      forces[2*j+1] -= Force[1];
	    }
	    j=Neighbours[2*j+1];
	  }
	}
	i=Neighbours[2*i+1];
      }
    }
  }
  free(Force);
}

#ifdef PRESSURE
void Move_Particles_ABP(particle* Particles,struct param Param,double* forces,double _time,long** Box,long* Neighbours,box *** NeighbouringBoxes, FILE* outputparam,double* pressure,double* average_density){
#else
void Move_Particles_ABP(particle* Particles,struct param Param,double* forces,double _time,long** Box,long* Neighbours,box *** NeighbouringBoxes, FILE* outputparam){
#endif
  long i;
  double newbi, newbj; // New boxes where particles will be moved
  double dr2,dr2_max;  // Used to compute the largest displacement
  double dx,dy;        // Used to compute displacements
  long i_max;          // index of particle with largest displacement
  
#ifdef LJ
  /*Compute force for each particle*/
  Force_LJ(Particles,Param,forces,Box,Neighbours,NeighbouringBoxes);
#endif

#ifdef WCA
  /*Compute force for each particle*/
  Force_LJ(Particles,Param,forces,Box,Neighbours,NeighbouringBoxes);
#endif

  /*Largest displacement*/
  dr2_max = 0;
  
  for(i=0;i<Param.N;i++){
    //Compute displacements along x and y
    
    dx = Param.v0 * cos(Particles[i].theta) * Param.dt + Param.sqrt2Dtdt*gasdevMT();
    dy = Param.v0 * sin(Particles[i].theta) * Param.dt + Param.sqrt2Dtdt*gasdevMT();

#ifdef LJ
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
    Particles[i].theta += Param.sqrt2Drdt*gasdevMT();
    /*
    if(Particles[i].theta>2*M_PI)
      Particles[i].theta-=2*M_PI;
    if(Particles[i].theta<0)
    Particles[i].theta+=2*M_PI;*/
    
    //Correct displacement for periodic boundary conditions
#ifdef PBC
    if(Particles[i].x>Param.Lx){
      Particles[i].x -= Param.Lx;
      if(Particles[i].x>Param.Lx){
	printf("x>2Lx\n");
	exit(1);
      }
    }
    if(Particles[i].x<0){
      Particles[i].x += Param.Lx;
      if(Particles[i].x<0){
	printf("x<-Lx\n");
	exit(1);
      }
    }
#endif
    //If closed, check you are not too close to the PBC
#ifdef CLOSEDBC
    if(Particles[i].x>Param.Lx-Param.sigma){
      printf("x_wall_right too large\n");
      exit(1);
    }
    if(Particles[i].x<Param.sigma){
      printf("x_wall_left too small\n");
      exit(1);
    }
#endif
    //Always use PBC along y
    if(Particles[i].y>Param.Ly){
      Particles[i].y -= Param.Ly;
      if(Particles[i].y>Param.Ly){
	printf("y>2Ly\n");
	exit(1);
      }
    }
    if(Particles[i].y<0){
      Particles[i].y += Param.Ly;
      if(Particles[i].y<0){
	printf("y<-Ly\n");
	exit(1);
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
  //If the largest displacement is larger than sigma/10, say so
  if(dr2_max>Param.sigma2*0.01){
    fprintf(outputparam,"at time %lg, dr_max=%lg largest acceptable:%lg\n",_time,sqrt(dr2_max),Param.sigma*0.1);
  }
}
 
double Distance2(particle* Particles,long i,long j,struct param Param){
  double dx=fabs(Particles[i].x-Particles[j].x);
  if(dx>Param.Lx/2.)
    dx=Param.Lx-dx;
  double dy=fabs(Particles[i].y-Particles[j].y);
  if(dy>Param.Ly/2.)
    dy=Param.Ly-dy;
  return dx*dx+dy*dy;
}
 
double SmallestDistance(particle* Particles,struct param Param){
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
