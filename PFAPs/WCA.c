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
#ifdef STRESSTENSOR
	   ,int StressBool,double** sigma_IK,double** sigma_a,long* count_measure_stress_tensor
#endif
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

#ifdef STRESSTENSOR
	    sigma_IK[bi][bj][0] +=dx*Force[0];
	    sigma_IK[bi][bj][1] +=dx*Force[1];
	    sigma_IK[bi][bj][2] +=dy*Force[0];
	    sigma_IK[bi][bj][3] +=dy*Force[1];
#endif
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
#ifdef STRESSTENSOR
	      //Box on top
	      if(k==1){
		//Coordinate of the intersection bet r_{ij} and the top of the box

		// r_{ij=(dx,dy)
		//  \tilde r_{ij}=(dx^*,dy^*)
		//   Colinearity tells us that
		//   dx/dy=dx^*/dy^* so that
		//  dx^*=dy^* dx/dy
		double yint=(bj+1)*Param.rbox;
		double dy1 = yint-Particles[i].y;
		double dx1 = dy1*dx/dy;
		double yint_partj=(nbj)*Param.rbox;
		double dy2 = Particles[j].y-yint_partj;
		double dx2 = dy2*dx/dy;
		
		sigma_IK[bi][bj][0] +=dx1*Force[0];
		sigma_IK[bi][bj][1] +=dx1*Force[1];
		sigma_IK[bi][bj][2] +=dy1*Force[0];
		sigma_IK[bi][bj][3] +=dy1*Force[1];
		sigma_IK[nbi][nbj][0] +=dx2*Force[0];
		sigma_IK[nbi][nbj][1] +=dx2*Force[1];
		sigma_IK[nbi][nbj][2] +=dy2*Force[0];
		sigma_IK[nbi][nbj][3] +=dy2*Force[1];
	      }
	      //Box on top right
	      if(k==2){
		// xint and yint are the right and top positions of particle-i's box
		double xint=(bi+1)*Param.rbox;
		double yint=(bj+1)*Param.rbox;
		// xint and yint are the bottom and left positions of particle-j's box
		double xint_partj=(nbi)*Param.rbox;
		double yint_partj=(nbj)*Param.rbox;
		
		double dx1 = xint-Particles[i].x;
		double dy1 = yint-Particles[i].y;		
		double dx2 = xint_partj-Particles[j].x;
		double dy2 = yint_partj-Particles[j].y;
		// abscissa of intersection of r_{ij} with tox of box i
		double dxa=dy1*dx/dy;		
		if(dxa>dx1){
		  dyp1=dx1*dy/dx;
		  dxa2=dxa-dx1;
		  dya2=dy1-dyp1;
		  dx2p=dx2-dxa2;
		  
		  sigma_IK[bi][bj][0] +=dx1*Force[0];
		  sigma_IK[bi][bj][1] +=dx1*Force[1];
		  sigma_IK[bi][bj][2] +=dyp1*Force[0];
		  sigma_IK[bi][bj][3] +=dyp1*Force[1];
		  
		  sigma_IK[nbi][nbj][0] +=dx2p*Force[0];
		  sigma_IK[nbi][nbj][1] +=dx2p*Force[1];
		  sigma_IK[nbi][nbj][2] +=dy2*Force[0];
		  sigma_IK[nbi][nbj][3] +=dy2*Force[1];

		  nnbi=NeighbouringBoxes[bi][bj][3].i;
		  nnbj=NeighbouringBoxes[bi][bj][3].j;
		  
		  sigma_IK[nnbi][nnbj][0] +=dxa2*Force[0];
		  sigma_IK[nnbi][nnbj][1] +=dxa2*Force[1];
		  sigma_IK[nnbi][nnbj][2] +=dya2*Force[0];
		  sigma_IK[nnbi][nnbj][3] +=dya2*Force[1];
		}
		else{
		  dxp1=dxa;
		  dy2p=dx2*dy/dx;
		  dya2=dy2-dy2p;
		  dxa2=dx1-dxa;
		  
		  sigma_IK[bi][bj][0] +=dxp1*Force[0];
		  sigma_IK[bi][bj][1] +=dxp1*Force[1];
		  sigma_IK[bi][bj][2] +=dy1*Force[0];
		  sigma_IK[bi][bj][3] +=dy1*Force[1];
		  
		  sigma_IK[nbi][nbj][0] +=dx2*Force[0];
		  sigma_IK[nbi][nbj][1] +=dx2*Force[1];
		  sigma_IK[nbi][nbj][2] +=dy2p*Force[0];
		  sigma_IK[nbi][nbj][3] +=dy2p*Force[1];

		  nnbi=NeighbouringBoxes[bi][bj][1].i;
		  nnbj=NeighbouringBoxes[bi][bj][1].j;
		  
		  sigma_IK[nnbi][nnbj][0] +=dxa2*Force[0];
		  sigma_IK[nnbi][nnbj][1] +=dxa2*Force[1];
		  sigma_IK[nnbi][nnbj][2] +=dya2*Force[0];
		  sigma_IK[nnbi][nnbj][3] +=dya2*Force[1];

		}
	      }
	      //Box on right
	      if(k==3){
		// r_{ij=(dx,dy)
		//  \tilde r_{ij}=(dx^*,dy^*)
		//   Colinearity tells us that
		//   dx/dy=dx^*/dy^* so that
		//  dy^*=dx^* dy/dx
		// xint is the right of the box of particle i
		double xint=(bi+1)*Param.rbox;
		double dx1 = xint-Particles[i].x;
		double dy1 = dx1*dy/dx;
		//xint_partj is the left of the box of particle j
		double xint_partj=(nbi)*Param.rbox;
		double dx2 = Particles[j].x-xint_partj;
		double dy2 = dx2*dy/dx;
		
		sigma_IK[bi][bj][0] +=dx1*Force[0];
		sigma_IK[bi][bj][1] +=dx1*Force[1];
		sigma_IK[bi][bj][2] +=dy1*Force[0];
		sigma_IK[bi][bj][3] +=dy1*Force[1];
		sigma_IK[nbi][nbj][0] +=dx2*Force[0];
		sigma_IK[nbi][nbj][1] +=dx2*Force[1];
		sigma_IK[nbi][nbj][2] +=dy2*Force[0];
		sigma_IK[nbi][nbj][3] +=dy2*Force[1];
	      }
	      //Box on bottom right: update
	      if(k==4){
		
	      }
#endif
	      
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
