#include "./bessel.c"
 
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




void GradKer_Bessel_Hashing_PBC(particle* Particles,struct param Param,double* GradKer,int** Box,int* Neighbours,box *** NeighbouringBoxes){
  int i,j,k,bi,bj,nbi,nbj;
  double dist,dist2;
  double amp;
  double dx,dy,epsilonx,epsilony;
 

  //Reset the array density at [0,0,0,...,0]
  memset(GradKer, 0, sizeof(double)*2*Param.N);
  
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
	  // if j contributes to GradKer[i], the converse is also true
	  // Here we use a TwoSchwartzShoulders
	  if(dist2<Param.rmax2){
	    dist=sqrt(dist2);

	    amp=-( Param.a1*Param.s1*bessk1(Param.s1*dist) + Param.a2*Param.s2*bessk1(Param.s2*dist) )/dist;
	    GradKer[2*i]   += amp*dx;
	    GradKer[2*i+1] += amp*dy;
	    GradKer[2*j]   -= amp*dx;// grad K(x_j-x_i) = - grad K(x_i-x_j)
	    GradKer[2*j+1] -= amp*dy;
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
	    // if j contributes to GradKer[i], the converse is also true
	    // Here we use a TwoSchwartzShoulders
	    if(dist2<Param.rmax2){ 
	      dist=sqrt(dist2);
	      
	      amp=-( Param.a1*Param.s1*bessk1(Param.s1*dist) + Param.a2*Param.s2*bessk1(Param.s2*dist) )/dist;
	      GradKer[2*i]   += amp*dx;
	      GradKer[2*i+1] += amp*dy;
	      GradKer[2*j]   -= amp*dx;// grad K(x_j-x_i) = - grad K(x_i-x_j)
	      GradKer[2*j+1] -= amp*dy;
	    }
	    j=Neighbours[2*j+1];
	  }
	}
	i=Neighbours[2*i+1];
      }
    }
  }
}





void Move_Particles_RTP(particle* Particles,struct param Param,double* GradKer,double* alpha,double _time,int** Box,int* Neighbours,box *** NeighbouringBoxes,double* dt_temp,double* NegTumbRateCount){
  int i;
  double ux, uy;
  double newbi, newbj;
  double vdt_temp, alpha_max;
  
  /*Compute GradKer for each particle*/
  GradKer_Bessel_Hashing_PBC(Particles,Param,GradKer,Box,Neighbours,NeighbouringBoxes);

  NegTumbRateCount[0] = 0; //At each time step we reset the number of neg. tumb. rates
  
  /* Compute the tumbling rate of each particle */
  for(i=0;i<Param.N;i++){
    
    //Compute orientation
    ux=cos(Particles[i].theta);
    uy=sin(Particles[i].theta);
    
    //update tumbling rate
    alpha[i] = Param.alpha0 + Param.alpha1 * ( ux*GradKer[2*i] + uy*GradKer[2*i+1] );

    //Put to zero any negative tumbling rate and had +1 to the counter of neg tumbling rate
    if(alpha[i]<=0)
      {
	alpha[i] = 0;
	NegTumbRateCount[0] +=1;
      }
  }


  
  /*Is it necessary to change dt_temp ?*/
  //Compute alpha_max
  alpha_max = 0;//initialize max tumbling rate to zero

  for(i=0;i<Param.N;i++){
    if (alpha[i] > alpha_max){
      alpha_max = alpha[i] ;
    }
  }
  
  //rescale the time-step in order to make the average number of tumble (of the most tumbling particle) alpha*dt_temp reasonable. Lets say of the order of 1/5.
  if (alpha_max* Param.dt > 1/5){
    dt_temp[0] = 1/( 5*alpha_max );
  }
  else{
    dt_temp[0]=Param.dt;
  }
  
  
  
  
  /* Compute change in orientation then displacement */
  for(i=0;i<Param.N;i++){

     //Test for tumbles
    if(alpha[i]>0) // if alpha[i]=0 then Particles[i].time tends to +infinity so no tumble
      {
	if(genrand64_real2() < alpha[i]*dt_temp[0]){
	  Particles[i].theta = 2.*M_PI*genrand64_real2();
	  //while(_time+dt_temp[0]>Particles[i].time){
	  //Particles[i].theta = 2.*M_PI*genrand64_real2();
	  //Particles[i].time  = Particles[i].time - 1./Param.alpha[i]*log(genrand64_real3());
	}
      }
    //Compute orientation again 
    ux=cos(Particles[i].theta);
    uy=sin(Particles[i].theta);
    
    //Move the particle
    vdt_temp =  Param.v0*dt_temp[0];
    Particles[i].x     += vdt_temp*ux;
    Particles[i].y     += vdt_temp*uy;

    //Correct displacement for periodic boundary conditions
#ifdef PBC
    Particles[i].x=(Particles[i].x<Param.Lx)?Particles[i].x:(Particles[i].x-Param.Lx);
    Particles[i].x=(Particles[i].x>=0)?Particles[i].x:(Particles[i].x+Param.Lx);
    Particles[i].y=(Particles[i].y<Param.Ly)?Particles[i].y:(Particles[i].y-Param.Ly);
    Particles[i].y=(Particles[i].y>=0)?Particles[i].y:(Particles[i].y+Param.Ly);
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
}

