
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>

int takestep_random_displacement(double * x, int n, double stepsize, gsl_rng * gslrng)
{
  int i = 0;
  stepsize *= 2;

  for (i=0; i<n; ++i){
    x[i] += stepsize * (gsl_rng_uniform(gslrng) - 0.5);
  }
}

int mc(double *x0, double *xreturn, int natoms, long int mciter, double stepsize, double Emax, double radius, long int seed)
{
  //set up the random number generator
  gsl_rng * gslrng; //declare the random number generator
  gsl_rng_env_setup(); //set up rng
  //const gsl_rng_type * T = gsl_rng_default; //declare type of random number generator
  const gsl_rng_type * T = gsl_rng_mt19937;
  gslrng = gsl_rng_alloc (T); //allocate working space for the rng
  gsl_rng_set (gslrng, seed);

  double eps = 1.;
  double sig = 1.;

  int N = natoms*3;
  int i, istep;
  int accept = 1;

  double *x = malloc(sizeof(double) * N);
  double *xnew = malloc(sizeof(double) * N);
  double *xtemp = NULL;

  double E, Enew;


  for (i=0; i<N; ++i){
    x[i] = x0[i];
  }

  E = ljenergy(x, N, eps, sig);


  printf("in mc.c\n");

  for (istep=0; istep<N; ++istep){

    for (i=0; i<N; ++i){
      xnew[i] = x[i];
    }
    
    takestep_random_displacement(xnew, 3*natoms, stepsize, gslrng);

    Enew = ljenergy(xnew, N, eps, sig);

    accept = ( Enew < Emax );
    if ( accept ){
      //test spherical container
      
      if (accept){
        xtemp = x;
        x = xnew;
        xnew = xtemp;
        E = Enew;
      }
    }



  }

  for (i=0; i<N; ++i){
    xreturn[i] = x[i];
  }

  free(x);
  free(xnew);
  gsl_rng_free (gslrng); //free the rng working space
}
