#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>



double label_energy(double E, double macheps, gsl_rng * r)
{
  double label, u;
  
  u = gsl_rng_uniform(r);
    
  label = 1 + (u - 0.5) * macheps;
  //printf("u %g\n", u);
  //printf("macheps %g\n", macheps);
  //printf("label %g\n", label);
  E = E * label;

  return E;
}

double get_energy_change(long int * neighbor_list, long int * nbegin, long int * nend, long int * spins, long int itrial)
{
  double E, deltaS;
  long int n0, n, nf;
  long int ineib;
  E = 0.;
  deltaS = - 2 * spins[itrial]; // = new_spin - old_spin

  n0 = nbegin[itrial];
  nf = nend[itrial];
  for (n = n0; n < nf; ++n)
    {
      ineib = neighbor_list[n];
      //printf("itrial ineib %d %d\n", itrial, ineib);
      E -= deltaS * spins[ineib];
    }
  return E;
}

long int mcising(long int *spins, double *Energy, long int nspins, long long int niter, double Emax, long long int seed,
		 long int *neighbor_list, long int *nbegin, long int *nend)
{
  long int istep;
  long int itrial;
  double deltaE, E, new_E;
  long int accept;
  long int naccept = 0;
  double macheps = 0.00000000001;
  //set up the random number generator
  gsl_rng * gslrng; //declare the random number generator
  gsl_rng_env_setup(); //set up rng
  //const gsl_rng_type * T = gsl_rng_default; //declare type of random number generator
  const gsl_rng_type * T = gsl_rng_mt19937;
  gslrng = gsl_rng_alloc (T); //allocate working space for the rng
  gsl_rng_set (gslrng, seed);

  E = *Energy;

  //printf("sizeof int, long int %d %d\n", sizeof(int), sizeof(long int));

  for (istep=0; istep<niter; ++istep)
    {
      itrial = gsl_rng_uniform_int(gslrng, nspins);
      
      deltaE = get_energy_change(neighbor_list, nbegin, nend, spins, itrial);
      // printf("deltaE %g\n", deltaE);
      
      new_E = E + deltaE;
      //printf("newE %g\n", new_E);
      accept = ( new_E <=  Emax ); //I have intentionally changed to strictly less
      
      if ( accept )
	{
	  spins[itrial] *= -1;
	  E = new_E;
	  naccept += 1;
	}
    }

  //*Energy = label_energy(E, macheps, gslrng);
  *Energy = E;
  gsl_rng_free (gslrng); //free the rng working space
  //printf("done compiled mc: %d %d %f stepsize %f n_energy_ok %d natoms %d\n", naccept, istep, E, stepsize, energy_ok_count, natoms);
  return naccept;
}
