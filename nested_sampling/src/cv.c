#include <math.h>
#include <stdlib.h>
#include <stdio.h>
double X_imp(int i, double K, double P);
void compute_dos(long double* gl, int N, double P, double K);
void renorm_energies(double* El, int N, double Emin);
double heat_capacity(double* El, long double* gl, int N, double T, double ndof);
void heat_capacity_loop(double* El, long double* gl, double* Cvl, int N, double Tmin, double Tmax, int nT, double ndof);

double X_imp(int i, double K, double P)
{
  double X;
  X = (K - ((i+1)%(int)P) )/( K - ((i+1)%(int)P) + 1);
  return X;
}

void compute_dos(long double* gl, int N, double P, double K)
{
  // gl is an array of 0's of size N, K is the number of replicas
  int i;
  long double X = 1 - P/(K+1);
  long double Xm = X; //this is X1
  long double Xf, Xb;
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

void compute_dos_imp(long double* gl, int N, double P, double K)
{
  // gl is an array of 0's of size N, K is the number of replicas
  int i;
  long double X;
  long double Xm = X_imp(0,K,P); //this is X1
  long double Xf, Xb;
  Xb = 2-X_imp(0,K,P); // reflecting boundary condition, this is X0
  Xf = Xm * X_imp(1,K,P); 
  gl[0] = 0.5 * (Xb - Xf);
  //calculate density of states for stored energies
  for(i=1;i<(N-K-1);++i)
    {
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
double heat_capacity(double* El, long double* gl, int N, double T, double ndof)
{
  //K is the number of replicas, beta the reduced temperature and E is the array of energies 
  int i;
  double Cv;
  long double bolz;
  long double Z = 0;
  long double U = 0;
  long double U2 = 0;
  long double beta = 1/T;

  for(i=0;i<N;++i)
  {
    bolz = gl[i]*exp(-beta*El[i]); 
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
void heat_capacity_loop(double* El, long double* gl, double* Cvl, int N, double Tmin, double Tmax, int nT, double ndof)
{
  //Cvl is a 0's array of size N (same size as El)
  int i;
  double dT = (Tmax - Tmin) / nT;
  double T = Tmin;
  
  for(i=0;i<nT;++i)
  {
    Cvl[i] = heat_capacity(El, gl, N, T, ndof);
    T += dT;
  }
}
