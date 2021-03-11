double Distance2_obst(obstacle* Obstacles,int i,int j,param Param){
  double dx=fabs(Obstacles[i].x-Obstacles[j].x);
  if(dx>Param.Lx/2.)
    dx=Param.Lx-dx;
  double dy=fabs(Obstacles[i].y-Obstacles[j].y);
  if(dy>Param.Ly/2.)
    dy=Param.Ly-dy;
  return dx*dx+dy*dy;
}

double Distance2_partobst(particle* Particles,long i,obstacle* Obstacles,int j,param Param){
  double dx=fabs(Particles[i].x-Obstacles[j].x);
  if(dx>Param.Lx/2.)
    dx=Param.Lx-dx;
  double dy=fabs(Particles[i].y-Obstacles[j].y);
  if(dy>Param.Ly/2.)
    dy=Param.Ly-dy;
  return dx*dx+dy*dy;
}

/*
  This function adds obstacle i to the Box[bi][bj]. It always add the
  obstacle as the first one and fix the links between obstacles.
*/
void AddinBox(int i, int bi,int bj,int** Box,obstacle* Obstacles){
  int k;
  
  if(Box[bi][bj]==-1){
    Box[bi][bj]=i;
    Obstacles[i].prev=-1;
    Obstacles[i].next=-1;
  }
  else{
    //The first obstacle was k
    k=Box[bi][bj];
    //The new first obstacle is i
    Box[bi][bj]=i;
    //The next obstacle following i is k
    Obstacles[i].next=k;
    //The obstacle before k is i
    Obstacles[k].prev=i;
    //There is no one before i
    Obstacles[i].prev=-1;
  }
}

void Place_Obstacles(obstacle* Obstacles, param Param, int** Box){
  int i,j;
  double mindist2,distance2;
  if(Param.Nobs==1){
    Obstacles[0].x           = genrand64_real3()*Param.Lx;
    Obstacles[0].y           = genrand64_real3()*Param.Ly;
    // Add obstacle in the good box.
    Obstacles[0].bi    = floor(Obstacles[0].x/Param.rbox);
    Obstacles[0].bj    = floor(Obstacles[0].y/Param.rbox);
    AddinBox(0,Obstacles[0].bi,Obstacles[0].bj,Box,Obstacles);
  }
  else{
    for(i=0;i<Param.Nobs;i++){
      // As long as the obstacles overlap, try a new position
      mindist2=0;
      while(mindist2<Param.sigma*Param.sigma){
	// Get a new position
	Obstacles[i].x           = genrand64_real3()*Param.Lx;
	Obstacles[i].y           = genrand64_real3()*Param.Ly;
	
	// Compute the smallest distance (squared) to existing obstacles
	mindist2=Param.Lx*Param.Ly;
	for(j=0;j<i;j++){
	  distance2 = Distance2_obst(Obstacles,i,j,Param);
	  if(distance2<mindist2)
	    mindist2=distance2;
	}
      }
      
      // Add obstacle in the good box.
      Obstacles[i].bi    = floor(Obstacles[i].x/Param.rbox);
      Obstacles[i].bj    = floor(Obstacles[i].y/Param.rbox);
      AddinBox(i,Obstacles[i].bi,Obstacles[i].bj,Box,Obstacles);
    }
  }
};

void Place_Particle(long i,particle* Particles, param Param,obstacle* Obstacles){
  
  // As long as the particle is in an obstacle, try a new position
  double mindist2,distance2;
  int j;
  
  mindist2=0;
  while(mindist2<Param.sigma*Param.sigma){
    // Get a new position
    
    Particles[i].x           = genrand64_real3()*Param.Lx;
    Particles[i].y           = genrand64_real3()*Param.Ly;

    mindist2=Param.Lx*Param.Ly;
    for(j=0;j<Param.Nobs;j++){
      distance2 = Distance2_partobst(Particles,i,Obstacles,j,Param);
      if(distance2<mindist2)
	mindist2=distance2;
    }
  }
  
  Particles[i].theta       = genrand64_real3()*2*M_PI;

  // This is a real number distributed according to \alpha exp(-\alpha t)
  Particles[i].next_tumble = -1 / Param.alpha * log(genrand64_real3()); 
  
  //Add particle in the good box.
  Particles[i].bi    = floor(Particles[i].x/Param.rbox);
  Particles[i].bj    = floor(Particles[i].y/Param.rbox);
};

/* 
   This function loops through all the neigbouring boxes of the box
   (bi,bj) and check for possible interactions between obstacles and a
   particle move to (x+dx) (y+dy). In the case of a collision,
   dr_\perp is set to zero.

   Question: are sequential collisions handled properly ?
*/

void Collisions_obstacles(double* dx,double* dy,double x,double y,obstacle* Obstacles,int bi,int bj,param Param,box *** NeighbouringBoxes,int** Box){
  int k; //counter for the box
  int _obstacle; //number of the obstacle
  int nbi,nbj; // indices of the box under scrutiny
  double epsilonx,epsilony,dxobs,dyobs,dist2,norm2,factor,xobs,yobs;
  
  for(k=0;k<9;k++){
    nbi       = NeighbouringBoxes[bi][bj][k].i;
    nbj       = NeighbouringBoxes[bi][bj][k].j;
    epsilonx  = NeighbouringBoxes[bi][bj][k].epsilonx;
    epsilony  = NeighbouringBoxes[bi][bj][k].epsilony;
    
    // Take the first obstacle of the neighbouring box
    _obstacle = Box[nbi][nbj];
    
    //Cycle through all particles after j
    while(_obstacle!=-1){
      xobs  = Obstacles[_obstacle].x + epsilonx;
      yobs  = Obstacles[_obstacle].y + epsilony;
      dxobs = fabs(xobs-x-dx[0]);
      dyobs = fabs(yobs-y-dy[0]);
      dist2 = dxobs*dxobs + dyobs*dyobs;
      // if the particles interact with the obstacle, set dr_perp to zero
      if(dist2<Param.sigma2){
	norm2=(xobs-x)*(xobs-x)+(yobs-y)*(yobs-y);
	factor = ( dx[0]*(xobs-x) + dy[0]*(yobs-y) )/norm2;
	dx[0] -= factor * (xobs-x);
	dy[0] -= factor * (yobs-y);	
      }
      _obstacle=Obstacles[_obstacle].next;
    }
  }
}

void Move_Particles_RTP(particle* Particles, param Param,double *_time,obstacle* Obstacles,int** Box,box*** NeighbouringBoxes){
  long i;
  int newbi,newbj;
  double dx,dy;
  _time[0] += Param.dt;
  
  for(i=0;i<Param.N;i++){
    dx = Param.v0 * cos( Particles[i].theta ) * Param.dt;
    dy = Param.v0 * sin( Particles[i].theta ) * Param.dt;
    
    Collisions_obstacles(&dx,&dy,Particles[i].x,Particles[i].y,Obstacles,Particles[i].bi,Particles[i].bj,Param,NeighbouringBoxes,Box);
    
    Particles[i].x += dx;
    Particles[i].y += dy;
    
    if(_time[0]>Particles[i].next_tumble){
      Particles[i].theta       = genrand64_real3()*2*M_PI;
      
      // This is a real number distributed according to \alpha exp(-\alpha t)
      Particles[i].next_tumble = _time[0] -1 / Param.alpha * log(genrand64_real3()); 
    }
    
#ifdef PERIODICBC
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
#endif
    //Test if the particle has changed box
    newbi    = floor(Particles[i].x/Param.rbox);
    newbj    = floor(Particles[i].y/Param.rbox);
    
    if(Particles[i].bi!=newbi||Particles[i].bj!=newbj){
      Particles[i].bi=newbi;
      Particles[i].bj=newbj;
    }
  }
};

void Record_positions(double _time, particle* Particles, param Param, FILE* outputpos){
  long i;
  // Store positions, angles, packing fraction using a particle diameter sigma
  for(i=0;i<Param.N;i++){
    fprintf(outputpos,"%lg\t%d\t%lg\t%lg\t%lg\n",_time,i,Particles[i].x,Particles[i].y,Particles[i].theta);
  }
  fprintf(outputpos,"\n");
  fflush(outputpos);
}

void Measure_density_current(particle* Particles,double** density,double*** current,param Param){
  int i;
  int bi,bj;
  
  for(i=0;i<Param.N;i++){
    bi=(int) floor(Particles[i].x/Param.dr);
    bj=(int) floor(Particles[i].y/Param.dr);
    density[bi][bj]    += 1;
    current[bi][bj][0] += Param.v0*cos(Particles[i].theta); // in principle, should be projected orthogonally to obstacles
    current[bi][bj][1] += Param.v0*sin(Particles[i].theta); // in principle, should be projected orthogonally to obstacles
  }
}

void Record_density_current(double** density,double*** current,FILE* outputdens,param Param,double densitycount,double _time){
  int i,j;
  
  for(i=0;i<Param.NxBoxDensity;i++){
    for(j=0;j<Param.NyBoxDensity;j++){
      fprintf(outputdens,"%lg\t%lg\t%lg\t",_time,(i+.5)*Param.dr,(j+.5)*Param.dr);
      fprintf(outputdens,"%lg\t",density[i][j]/densitycount/Param.dr/Param.dr);
      fprintf(outputdens,"%lg\t%lg\n",current[i][j][0]/densitycount/Param.dr/Param.dr,current[i][j][1]/densitycount/Param.dr/Param.dr);
      density[i][j]=0;
      current[i][j][0]=0;
      current[i][j][1]=0;
    }
    fprintf(outputdens,"\n");
  }
}

