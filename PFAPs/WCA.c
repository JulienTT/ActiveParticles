/*
  This function finds dt such that

  orderdt^2 dt^2 + orderrootdt^2 dt = rmax
  
 */
double FindDt(particle* Particles,double* orderdt,double* orderrootdt,long i_max,param Param){
  double a,b;
  a = orderdt[2*i_max]*orderdt[2*i_max] + orderdt[2*i_max+1]*orderdt[2*i_max+1];
  b = orderrootdt[2*i_max]*orderrootdt[2*i_max] + orderrootdt[2*i_max+1]*orderrootdt[2*i_max+1];
  return (-b + sqrt(b*b+4.*a*Param.rmax2) ) / (2.*a);
};

void WCA_Force(double* Force,double dx,double dy,double dist2,param Param){
  double dist6=dist2*dist2*dist2;
  double amp=Param.amp/(dist6*dist2)*(1-2*Param.sigma6/dist6);
  Force[0]=dx*amp;
  Force[1]=dy*amp;
}

double LJ_Potential(double dist2,param Param){
  double sigma6overdist6=Param.sigma6/(dist2*dist2*dist2);
  return 4*Param.epsilon*sigma6overdist6*(sigma6overdist6-1)+Param.epsilon;
}

void Force_LJ(particle* Particles,param Param, double* forces,long** Box,long* Neighbours,box*** NeighbouringBoxes
#ifdef OBSERVABLE
	      , double* Energy
#endif
	      ){
  
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
      // their contribution to the forces
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
#ifdef OBSERVABLE
	    //This compute the interaction energy between parts i and j
	    Energy[0] += LJ_Potential(dist2,Param)*Param.dt;
#endif
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
#ifdef OBSERVABLE
	      //This computes the interaction energy between parts i and j
	      Energy[0] += LJ_Potential(dist2,Param)*Param.dt;
#endif
	      //This computes the force on particle i
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
