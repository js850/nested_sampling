#include <math.h>
#include <stdio.h>

// TODO: make a lennard jones class! for now just do a function for testing

double lj(double *x1, double *x2, double eps, double sig)
{
  double r=0, r2=0;
  double sig2 = sig*sig;
  double ir6 = 0;
  int i;
  
  for(i=0; i<3; ++i) {
    r = x2[i] - x1[i];
    r2 += r*r;
  }
  ir6 = pow(sig2/r2, 3);
  
  return 4.*eps*(ir6 * ir6 - ir6);
}

double ljenergy(double *x, int N, double eps, double sig)
{
  int i, j;
  double energy = 0;
  for(i=0; i<N; i+=3) 
    for(j=i+3; j<N; j+=3)
      energy += lj(&x[i], &x[j], eps, sig);
  printf("energy %f ", energy);
  return 11.;
  return energy;
}

