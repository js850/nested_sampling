#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
double max_array(double* a, int N);
double X_imp(int i, double K, double P);
void compute_dos(double* gl, int N, double P, double K);
void renorm_energies(double* El, int N, double Emin);
void log_weigths(double* El, double* gl, double* wl, int N, double T);
double heat_capacity(double* El, double* gl, int N, double T, double ndof);
void heat_capacity_loop(double* El, double* gl, double* wl, double* Cvl, int N, double Tmin, double Tmax, int nT, double ndof);

double max_array(double* a, int N)
{
  int i;
  double max=-DBL_MAX;
  for (i=0; i<N; ++i)
    {
      if (a[i]>max)
	{
	  max=a[i];
	}
    }
  return max;
}

double X_imp(int i, double K, double P)
{
  double X;
  X = (K - ((i+1)%(int)P) )/( K - ((i+1)%(int)P) + 1);
  return X;
}

void compute_dos(double* gl, int N, double P, double K)
{
  // gl is an array of 0's of size N, K is the number of replicas
  int i;
  double X = 1 - P/(K+1);
  double Xm = X; //this is X1
  double Xf, Xb;
  Xb = 2-X; // reflecting boundary condition, this is X0
  Xf = Xm * X; 
  gl[0] = 0.5 * (Xb - Xf);
  for(i=1;i<(N-1);++i)
    {
      Xb = Xm;
      Xm = Xf;
      Xf = Xm * X;
      gl[i] = 0.5 * (Xb - Xf);
    }
  Xb = Xm;
  Xm = Xf;
  Xf = -Xm;
  gl[N-1] =  0.5 * (Xb - Xf);  
}

void compute_dos_imp(double* gl, int N, double P, double K)
{
  // gl is an array of 0's of size N, K is the number of replicas
  int i;
  double X;
  double Xm = X_imp(0,K,P); //this is X1
  double Xf, Xb;
  Xb = 2-X_imp(0,K,P); // reflecting boundary condition, this is X0
  Xf = Xm * X_imp(1,K,P); 
  gl[0] = 0.5 * (Xb - Xf);
  //calculate density of states for stored energies
  for(i=1;i<(N-K-1);++i)
    {
      //printf("Xm %E \n",Xm);
      Xb = Xm;
      Xm = Xf;
      Xf = Xm * X_imp(i+1,K,P);
      gl[i] = 0.5 * (Xb - Xf);
    }
  //calculate density of states for live replica energies
  for(i=(N-K-1);i<(N-1);++i)
    {
      Xb = Xm;
      Xm = Xf;
      Xf = Xm * (K-i)/(K-i+1);
      gl[i] = 0.5 * (Xb - Xf);
    }
  Xb = Xm;
  Xm = Xf;
  Xf = -Xm;
  gl[N-1] =  0.5 * (Xb - Xf);  
}

void log_weigths(double* El, double* gl, double* wl, int N, double T)
{
  int i;
  double beta = 1/T;
  for(i=0;i<N;++i)
    {
      wl[i] = log(gl[i]) - beta * El[i];
    }
}

////////////renormalise energies wrt ground state/////////////////
void renorm_energies(double* El, int N, double Emin)
{
  int i;
  for(i=0;i<N;++i)
  {
    El[i] -= Emin;
  }
}

////////////////////////////caclulate heat capacity for a single T////////////////////////
double heat_capacity(double* El, double* wl, int N, double T, double ndof)
{
  //K is the number of replicas, beta the reduced temperature and E is the array of energies 
  int i;
  double Cv;
  double bolz;
  double Z = 0;
  double U = 0;
  double U2 = 0;
  double beta = 1/T;

  for(i=0;i<N;++i)
  {
    bolz = exp(wl[i]); 
    Z += bolz;
    U += El[i] * bolz;
    U2 += El[i] * El[i] * bolz;
  }

  U /= Z;
  U2 /= Z;
  Cv =  (U2 - U*U)*beta*beta + ndof/2;
  
  return Cv;
}

//////////////////////////////calculate heat capacity over a set of Ts/////////////////////
void heat_capacity_loop(double* El, double* gl, double* wl, double* Cvl, int N, double Tmin, double Tmax, int nT, double ndof)
{
  //Cvl is a 0's array of size N (same size as El)
  int i,j;
  double dT = (Tmax - Tmin) / nT;
  double T = Tmin;
  double wl_max;
  
  for(i=0;i<nT;++i)
  {
    log_weigths(El, gl, wl, N, T);
    wl_max = max_array(wl,N);
    
    //printf("wl_max %d \n",wl_max);
    
    for(j=0;j<N;++j)
      {
	wl[j] -= wl_max;
      }
    
    Cvl[i] = heat_capacity(El, wl, N, T, ndof);
    T += dT;
  }
}
