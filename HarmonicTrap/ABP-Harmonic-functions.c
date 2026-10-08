void ERROR(char* msg){
  printf(msg);
  exit(1);
}

void Move_Particles_ABP(particle* Particles,param Param){
  long i;
  double dx,dy;
  
  for( i=0 ; i<Param.N ; i++){
    dx = ( Param.v0 * cos(Particles[i].theta)  - Param.k* Particles[i].x ) * Param.dt + Param.sqrt2Dtdt*gasdev();
    dy = ( Param.v0 * sin(Particles[i].theta)  - Param.k* Particles[i].y ) * Param.dt + Param.sqrt2Dtdt*gasdev();
    Particles[i].x += dx;
    Particles[i].y += dy;
    Particles[i].theta +=  Param.sqrt2Drdt*gasdev();
  }
}

/*Check input parameters, attribute them and store them*/
void  CheckAndReadInput(int argc, char* argv[], param* Param,
#ifdef _MT
			long long* seed,
#endif
#ifdef _PCG
			pcg128_t* seed,
#endif
			double* FinalTime, double* EquilibTime, histo* Histo, FILE** outputparam){
  int i;
  int argctarget=15;
  char command_base[1000]=""; // string that contains the desired format of the command line
  
  char name[200];      // string in which the file names are written
  strcat(command_base, "usage: ");
  strcat(command_base, argv[0]);
  strcat(command_base," file seed N dt v0 Dr Dt k FinalTime EquilibTime StoreHistoInter UpdateHistoInter dr rmax");
  
  if(argc!=argctarget){
    printf("%s\n",command_base);
    exit(1);
  }
  
  i=1;
  //File where parameters are stored
  sprintf(name,"%s-param",argv[i]);
  outputparam[0]=fopen(name,"w");
  
  //File where histogram is stored
  sprintf(name,"%s-histo",argv[i]);
  Histo[0].outputhisto=fopen(name,"w");
  
  i++;
#ifdef _MT
  seed[0]              = (long long) strtod(argv[i], NULL); i++;
#endif
#ifdef _PCG
  seed[0]              = (pcg128_t) strtod(argv[i], NULL); i++;
#endif
  
  Param[0].N           = (long) strtod(argv[i], NULL); i++;  
  Param[0].dt          = strtod(argv[i], NULL); i++;
  Param[0].v0          = strtod(argv[i], NULL); i++;
  Param[0].Dr          = strtod(argv[i], NULL); i++;
  Param[0].Dt          = strtod(argv[i], NULL); i++;
  Param[0].k           = strtod(argv[i], NULL); i++; 
  
  Param[0].sqrt2Drdt   = sqrt(2*Param[0].Dr*Param[0].dt);
  Param[0].sqrt2Dtdt   = sqrt(2*Param[0].Dt*Param[0].dt);  
  Param[0].Da          = Param[0].v0*Param[0].v0 / (2 * Param[0].Dr);
  
  FinalTime[0]         = strtod(argv[i], NULL); i++;
  EquilibTime[0]       = strtod(argv[i], NULL); i++;
  
  Histo[0].StoreHistoInter   = strtod(argv[i], NULL); i++;
  Histo[0].UpdateHistoInter  = strtod(argv[i], NULL); i++;  
  Histo[0].dr                = strtod(argv[i], NULL); i++;
  Histo[0].rmax              = strtod(argv[i], NULL); i++;
  Histo[0].NextStoreHisto    = Histo[0].StoreHistoInter;
  Histo[0].NextUpdateHisto   = Histo[0].UpdateHistoInter;  
  Histo[0].Nbin              = (int) (Histo[0].rmax/Histo[0].dr);
  Histo[0].histogram         = (double*) calloc(Histo[0].Nbin,sizeof(double));
  Histo[0].histocount        = 0;  
  
  /* Store parameters */
  fprintf(outputparam[0],"%s\n",command_base);
  
  for(i=0;i<argc;i++){
    fprintf(outputparam[0],"%s ",argv[i]);
  }
  fprintf(outputparam[0],"\n");
  
  printf("Max value for long: %ld\n", LONG_MAX);
  printf("Max value for int: %d\n", INT_MAX);
  
#ifdef _MT	     
  fprintf(outputparam[0],"seed is %lld\n", seed[0]);
#endif		     
#ifdef _PCG	     
  fprintf(outputparam[0],"seed is %" PRIu64 "\n", (uint64_t) seed[0]);
#endif		     

  fprintf(outputparam[0],"N is %ld\n", Param[0].N);
  fprintf(outputparam[0],"dt is %lg\n", Param[0].dt);
  fprintf(outputparam[0],"v0 is %lg\n", Param[0].v0);
  fprintf(outputparam[0],"Dr is %lg\n", Param[0].Dr);
  fprintf(outputparam[0],"Dt is %lg\n", Param[0].Dt);
  fprintf(outputparam[0],"k is %lg\n", Param[0].k);
  fprintf(outputparam[0],"Da is %lg\n", Param[0].Da);     
  
  fprintf(outputparam[0],"FinalTime is %lg\n",FinalTime[0]);
  fprintf(outputparam[0],"EquilibTime is %lg\n", EquilibTime[0]);
  fprintf(outputparam[0],"StoreInterHisto is %lg\n", Histo[0].StoreHistoInter);
  fprintf(outputparam[0],"UpdateInterHisto is %lg\n", Histo[0].UpdateHistoInter);  
  fprintf(outputparam[0],"binwidth (dr) is %lg\n", Histo[0].dr);
  fprintf(outputparam[0],"rmax is %lg\n", Histo[0].rmax);
  fprintf(outputparam[0],"Bin number is %ld\n", Histo[0].Nbin);  
  fflush(outputparam[0]);
}

void UpdateHisto(param Param, histo* Histo, particle* Particles){
  long i;
  double r;
  
  for(i=0 ; i<Param.N ; i++){
    r=sqrt( Particles[i].x * Particles[i].x + Particles[i].y * Particles[i].y );
    if (r<Histo[0].rmax){
      Histo[0].histogram[(int) (r/Histo[0].dr)] += 1.;
    }
    Histo[0].histocount += 1.0;
  }
  Histo[0].NextUpdateHisto += Histo[0].UpdateHistoInter;
}

void StoreHisto(param Param, histo* Histo,double _time){
  long i;
  for( i=0 ; i<Histo[0].Nbin ; i++ ){
    fprintf(Histo[0].outputhisto,"%lg\t%lg\t%lg\n",_time,( i + .5 ) * Histo[0].dr, \
	    Histo[0].histogram[i]/( Histo[0].histocount * M_PI * Histo[0].dr * Histo[0].dr * (2* i + 1) ) );
  }
  Histo[0].NextStoreHisto += Histo[0].StoreHistoInter;
  fprintf(Histo[0].outputhisto,"\n");
  fflush(Histo[0].outputhisto);
}
