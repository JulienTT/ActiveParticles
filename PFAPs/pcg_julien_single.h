/* generates a random number on [0,1]-real-interval */
double genrand64s_pcg_real1(pcg64s_random_t* rng)
{
  return (pcg64s_random_r(rng)>> 11) / 9007199254740991.0;
}

/* generates a random number on [0,1)-real-interval */
double genrand64s_pcg_real2(pcg64s_random_t* rng)
{
  return (pcg64s_random_r(rng) >> 11) / 9007199254740992.0;
}

/* generates a random number on (0,1)-real-interval */
double genrand64s_pcg_real3(pcg64s_random_t* rng)
{
  return ((pcg64s_random_r(rng) >> 12) + 0.5) / 4503599627370496.0;
  //return ldexp(((pcg64_random_r(rng) >> 12) + 0.5),-52);
  //return 0x1.0p-52 * ((pcg64_random_r(rng) >> 12) + 0.5);
}

int MEMORY=0;
double vMEM;

double gasdevPCG(pcg64s_random_t* rng){
  double fac,rsq,v1,v2;
  
  if (MEMORY == 0) {
    do {
      v1=2.0*genrand64s_pcg_real2(rng)-1.0;
      v2=2.0*genrand64s_pcg_real2(rng)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    
    vMEM=v1*fac;
    MEMORY=1;
    return v2*fac;
  } else {
    MEMORY=0;
    return vMEM;
  }
}

void gasdevPCG2(pcg64s_random_t* rng,double* Gauss){
  static double fac,rsq,v1,v2;
  
  do {
    v1=2.0*genrand64s_pcg_real2(rng)-1.0;
    v2=2.0*genrand64s_pcg_real2(rng)-1.0;
    rsq=v1*v1+v2*v2;
  } while (rsq >= 1.0 || rsq == 0.0);
  fac=sqrt(-2.0*log(rsq)/rsq);
  
  Gauss[0]=v1*fac;
  Gauss[1]=v2*fac;
}

