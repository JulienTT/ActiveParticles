#include "./bessel.c"

/*
  To add someone in a box, you need to put it there and to set its new
  neighbours.
*/
void AddinBoxOld(int i,int bi, int bj,int** Box,int* Neighbours){
  int k;
  //if there was nobody in the box, just add particle i
  if(Box[bi][bj]==-1){
    Box[bi][bj]=i;
    Neighbours[2*i]=-1;
    Neighbours[2*i+1]=-1;
  }
  // Otherwise add it at the end
  else{
    k=Box[bi][bj];
    while(Neighbours[2*k+1]!=-1){
      k=Neighbours[2*k+1];
    }
    Neighbours[2*k+1]=i;
    Neighbours[2*i]=k;
    Neighbours[2*i+1]=-1;
  }
}

/*
  To add someone in a box, you need to put it there and to set its new
  neighbours. This functions add at the top.
*/
void AddinBox(int i,int bi, int bj,int** Box,int* Neighbours){
  int k;
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
void RemovefromBox(int i,int bi,int bj,int** Box,int* Neighbours){
  int next;
  int prev;
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

void Density_TopHat_ClosedBC(particle* Particles,struct param Param,double* density){
  int i,j;
  double dist2;
  //Reset the array density at [0,0,0,...,0]
  memset(density, 0, sizeof(double)*Param.N);
  
  for(i=0;i<Param.N;i++){
    // reset the density
    // cycle through the particles to compute their contribution
    for(j=i+1;j<Param.N;j++){
      dist2=(Particles[i].x-Particles[j].x)*(Particles[i].x-Particles[j].x)+(Particles[i].y-Particles[j].y)*(Particles[i].y-Particles[j].y);
      // if j contributes to density[i], the converse is also true
      if(dist2<Param.rmax2){
	density[i] +=1.;
	density[j] +=1.;
      }
    }
    density[i] = density[i]/(M_PI*Param.rmax2);
  }
}

void Density_TopHat_PBC(particle* Particles,struct param Param,double* density){
  int i,j;
  double dist2;
  double dx,dy;
  //Reset the array density at [0,0,0,...,0]
  memset(density, 0, sizeof(double)*Param.N);
  
  for(i=0;i<Param.N;i++){
    // reset the density
    // cycle through the particles to compute their contribution
    for(j=i+1;j<Param.N;j++){
      // Check if particle is closer to image with PBC rather than
      // real particle
      dx=fabs(Particles[i].x-Particles[j].x);
      dy=fabs(Particles[i].y-Particles[j].y);
      dx=(dx<.5*Param.Lx)?dx:(Param.Lx-dx); // (test)?(if test true):(if test false)
      dy=(dy<.5*Param.Ly)?dy:(Param.Ly-dy);
      
      dist2=dx*dx+dy*dy;
      // if j contributes to density[i], the converse is also true
      // Here we use a top-hat kernel
      if(dist2<Param.rmax2){
	density[i] += 1.;
	density[j] += 1.;
      }
    }
    density[i] = density[i]/(M_PI*Param.rmax2);
  }
}

void Density_TopHat_Hashing_PBC(particle* Particles,struct param Param,double* density,int** Box,int* Neighbours,box *** NeighbouringBoxes){
  int i,j,k,bi,bj,nbi,nbj;
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
	  if(dist2<Param.rmax2){
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
	    if(dist2<Param.rmax2){
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
    density[i] = density[i]/(M_PI*Param.rmax2);
}








void Force_Bessel_Hashing_PBC(particle* Particles,struct param Param,double* force,int** Box,int* Neighbours,box *** NeighbouringBoxes){
  int i,j,k,bi,bj,nbi,nbj;
  double dist,dist2;
  double amp;
  double dx,dy,epsilonx,epsilony;
  

 

  //Reset the array density at [0,0,0,...,0]
  memset(force, 0, sizeof(double)*2*Param.N);
  
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
	  dx        = Particles[i].x-Particles[j].x;
	  dy        = Particles[i].y-Particles[j].y;
	  dist2     = dx*dx+dy*dy;
	  
	  /* if j contributes to force[i], the converse is also true. Here we use superposition of (linear combination of) green functions of screened Poisson eq.*/
	  if(dist2<Param.rmax2){
	    dist=sqrt(dist2);
	    
	     amp=-( Param.a1*Param.s1*bessk1(Param.s1*dist) + Param.a2*Param.s2*bessk1(Param.s2*dist) )/dist;
	     //GradV=(amp*dx,amp*dy) but force=-GradV
	    
	    force[2*i]   -= amp*dx;
	    force[2*i+1] -= amp*dy;
	    force[2*j]   += amp*dx;//Action-reaction principle
	    force[2*j+1] += amp*dy;
	    
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
	    dx=Particles[i].x-epsilonx-Particles[j].x;
	    dy=Particles[i].y-epsilony-Particles[j].y;
	    dist2=dx*dx+dy*dy;
	    // if j contributes to force[i], the converse is also true
	    
	    if(dist2<Param.rmax2){
	      dist=sqrt(dist2);

	      amp=-( Param.a1*Param.s1*bessk1(Param.s1*dist) + Param.a2*Param.s2*bessk1(Param.s2*dist) )/dist;
	      
	      force[2*i]   -= amp*dx;
	      force[2*i+1] -= amp*dy;
	      force[2*j]   += amp*dx;
	      force[2*j+1] += amp*dy;
	    }
	    j=Neighbours[2*j+1];
	  }
	}
	i=Neighbours[2*i+1];
      }
    }
  }
}


void Move_Particles_RTP(particle* Particles,struct param Param,double* force,double _time,int** Box,int* Neighbours,box*** NeighbouringBoxes,double* dt_temp){
  int i,newbi,newbj;
  double sqrt2Tdt_temp;//noise amplitude adapted to the new time step;
  double force_2_max;
  double force_max;

  
  /*Compute density via kernel*/
 Force_Bessel_Hashing_PBC(Particles,Param,force,Box,Neighbours,NeighbouringBoxes);

  
  /*Compute the maximum force amplitude and adapt the time step if necessary*/
  //initialize max force to zero
  force_2_max = 0;
  force_max   = 0;

    for(i=0;i<Param.N;i++)
    {
      //find the maximum force amplitude
      if(force_2_max < (force[2*i]*force[2*i] + force[2*i+1]*force[2*i+1]) )
	{
	  force_2_max = force[2*i]*force[2*i] + force[2*i+1]*force[2*i+1];
	}
    }
  
    //rescale the time-step in order to make the biggest (deterministic) move smaller than rmax/5
    force_max   = sqrt(force_2_max);
    if (force_max*Param.dt > Param.rmin/5)
      {
	dt_temp[0] = Param.rmin/(5*force_max);
      }
    else
      {
	dt_temp[0] = Param.dt;
      }
  
    if(dt_temp[0]!=Param.dt)
      {
	sqrt2Tdt_temp = sqrt(2.*dt_temp[0]/Param.beta);
      }
    else
      {
	sqrt2Tdt_temp = Param.sqrt2Tdt;
	//printf("pas de modif de dt time %lg\n",_time);
      }



  
  /* Compute displacement */
  for(i=0;i<Param.N;i++){
    //Move the particle
    Particles[i].x     += dt_temp[0]*force[2*i]+ sqrt2Tdt_temp*gasdevMT();
    Particles[i].y     += dt_temp[0]*force[2*i+1]+ sqrt2Tdt_temp*gasdevMT();
    
    //Correct displacement for periodic boundary conditions
#ifdef PBC
    Particles[i].x=(Particles[i].x<Param.Lx)?Particles[i].x:(Particles[i].x-Param.Lx);
    Particles[i].x=(Particles[i].x>=0)?Particles[i].x:(Particles[i].x+Param.Lx);
    Particles[i].y=(Particles[i].y<Param.Ly)?Particles[i].y:(Particles[i].y-Param.Ly);
    Particles[i].y=(Particles[i].y>=0)?Particles[i].y:(Particles[i].y+Param.Ly);

    //Test if the particle has changed box
    newbi    = floor(Particles[i].x/Param.rbox);
    newbj    = floor(Particles[i].y/Param.rbox);

    if(Particles[i].bi!=newbi||Particles[i].bj!=newbj){
      RemovefromBox(i,Particles[i].bi,Particles[i].bj,Box,Neighbours);
      AddinBox(i,newbi,newbj,Box,Neighbours);
      Particles[i].bi=newbi;
      Particles[i].bj=newbj;
    }
#endif
  }
}

