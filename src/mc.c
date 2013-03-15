#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>

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
  //printf("energy %f ", energy);
  //return 11.;
  return energy;
}

int takestep_random_displacement(double * x, int N, double stepsize, gsl_rng * gslrng)
{
  int i = 0;
  stepsize *= 2;

  for (i=0; i<N; ++i){
    x[i] += stepsize * (gsl_rng_uniform(gslrng) - 0.5);
  }
}

int check_spherical_container(double * x, int N, double radius)
{
  double radius2 = radius * radius;
  double r2;
  int natoms = N/3;
  int i, j;

  for (i=0; i<natoms; ++i){
    r2 = 0;
    for (j=0; j<3; ++j){
      r2 += x[3*i+j] * x[3*i+j];
    }
    if (r2 >= radius2){
      return 0;
    }
  }
  return 1;
}

int mc(double *x0, double *xreturn, int natoms, long int niter, double stepsize, double Emax, double radius, long int seed)
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
  int i;
  long int istep;
  int accept = 1;
  int naccept = 0;

  double *x = (double *) malloc(sizeof(double) * N);
  double *xnew = (double *) malloc(sizeof(double) * N);
  double *xtemp = NULL;

  double E, Enew;


  for (i=0; i<N; ++i){
    x[i] = x0[i];
  }

  E = ljenergy(x, N, eps, sig);
  //printf(" E %f Emax %f\n", E, Emax);


  //printf("in mc.c\n");

  for (istep=0; istep<niter; ++istep){

    for (i=0; i<N; ++i){
      xnew[i] = x[i];
    }
    
    takestep_random_displacement(xnew, N, stepsize, gslrng);

    Enew = ljenergy(xnew, N, eps, sig);
    //printf(" Enew %f Emax %f\n", Enew, Emax);
    //break;

    accept = ( Enew < Emax );
    if ( accept ){
      //test spherical container
      accept = check_spherical_container(x, N, radius);
      
      if (accept){
        for (i=0; i<N; ++i){
          x[i] = xnew[i];
        }
        //xtemp = x;
        //x = xnew;
        //xnew = xtemp;
        E = Enew;
        naccept += 1;
      }
    }

  }

  for (i=0; i<N; ++i){
    xreturn[i] = x[i];
  }

  free(x);
  free(xnew);
  gsl_rng_free (gslrng); //free the rng working space
  printf("%d %d %f\n", naccept, istep, E);
  return naccept;
}
