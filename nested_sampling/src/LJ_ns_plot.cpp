/* LJ_ns_plot.cpp
 *
 * Copyright (C) 2013 Stefano Martiniani <stefano.martiniani@gmail.com> 
 * at the University of Cambridge
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <iostream>
#include <chrono>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <getopt.h>

//#include "gsl/gsl_sf_gamma.h"

int main(int argc, char ** argv)
{
  ////////////////////////////////////OpenOutput/////////////////////////
  
  //const std::string file_name = "../../NS_gmin/lj31.energies";
  
  //Parser
  
  double T_init = 0.01;
  double T_fin = 1;
  double T_step = 0.000002;
  int opt_char;
  int K = 0;
  int L = 0;
  int P = 0;
  int n = 0;
  char* input_str = NULL;
  char* output_str = NULL;
  
  while ((opt_char = getopt (argc, argv, "I:O:L:K:n:P:i:f:s:h")) != EOF)
    switch (opt_char)
      {
      case 'I':
	input_str = optarg;
	break;
      case 'O':
	output_str = optarg;
	break;
      case 'L':
	L = atoi(optarg);
        break;
      case 'K':
	K = atoi(optarg);
        break;
      case 'n':
	n = atoi(optarg);
        break;
      case 'P':
	P = atoi(optarg);
        break;
      case 'i':
	T_init = atof(optarg);
	break;
      case 'f':
	T_fin = atof(optarg);
	break;
      case 's':
	T_step = atof(optarg);
	break;
      case 'h':
	std::cout<<"Options:"<<std::endl;
	std::cout<<"-I \t input file e.g /path/to/input.dat"<<std::endl;
	std::cout<<"-O \t output directory e.g. /path/to/directory "<<std::endl;
	std::cout<<"-L \t size of LJ cluster"<<std::endl;
	std::cout<<"-K \t number of replicas"<<std::endl;
	std::cout<<"-n \t number of mc steps (only needed to name directory) or number of degrees of freedom for L=1"<<std::endl;
	std::cout<<"-P \t number of processors"<<std::endl;
	std::cout<<"-i \t temperature lower bound, default: 0.02"<<std::endl;
	std::cout<<"-f \t temperature upper bound, default: 1"<<std::endl;
	std::cout<<"-s \t temperature step, default: 2e-5"<<std::endl;
	std::cout<<"-h \t help"<<std::endl;
	return 0;
      case '?':
        if (optopt == 'opt_char')
          fprintf (stderr, "Option -%opt_char requires an argument, use -h option for more info.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%opt_char', use -h option for more info.\n", optopt);
        else
          fprintf (stderr,"Unknown option character `-%opt_char', use -h option for more info.\n",optopt);
        return 1;
      default:
	abort ();
      }
  
  ///////////////////////////////////////////////////////////////                                                          
  std::string K_str = std::to_string(K);
  std::string L_str = std::to_string(L);
  std::string n_str = std::to_string(n);
  std::string directory = output_str;

  directory += "/LJ" + L_str + "/";
  mkdir(directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  directory += "K" + K_str + "_n" + n_str + "/";
  std::cout<<directory<<std::endl;
  mkdir(directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
  std::ofstream output;
  const std::string dos = directory + "dos.dat";
  const std::string free_energy = directory + "free_energy.dat";
  const std::string internal_energy = directory + "internal_energy.dat";
  const std::string heat_capacity = directory + "heat_capacity.dat";
  const std::string entropy = directory + "entropy.dat";
  const std::string sampled_shrinkage = directory + "sampled_shrinkage.dat";
  const std::string geometric_shrinkage = directory + "geometric_shrinkage.dat";
  const std::string arithmetic_shrinkage = directory + "arithmetic_shrinkage.dat";
  const std::string canonical_distribution2 = directory + "canonical_distribution1.dat";
  const std::string canonical_distribution22 = directory + "canonical_distribution95.dat";
  const std::string canonical_distribution24 = directory + "canonical_distribution15.dat";
  
  ///////////////////////////////////////////////////////////////////////
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
  /////////////////////////////SetConstants/////////////////////////////
  
  std::vector<double> En, Ena, list_X;
  list_X.reserve(K);
    
  /////////////////////////Read-datafile//////////////////////////////////
  int i;
  double T;
  std::ifstream chain_file;
  chain_file.open(input_str);

  if (chain_file.is_open())
    {
      while(chain_file>>T)
        {
	  En.push_back(T);
	  ++i;
        }
      chain_file.close();
    }  
  
  std::vector<double>::const_iterator itp;

  /*
    for (itp = En.begin(); itp != En.end(); itp++)
    {
      std::cout<< *itp << " ";
    }
  std::cout<<""<<std::endl;
  
  */
  
  //create list_X
  
  double const Kd = (double) K;
  int imax = En.size();
  double alpha = 1;

  i = 0;
  int lmax = En.size()-1;
  std::cout<<"size:"<<imax<<std::endl;
  double min, max, ni, mi;
  
  if (L == 0)
    {
      mi = 0;
      while (i <= lmax)
	{
	  ni = 1;
	  Ena.push_back(En.at(i));
	  
	  if ((i+1) < lmax)
	    {
	      min = abs(En.at(i+1)) - 0.1;
	      max = abs(En.at(i+1)) + 0.1;
	      
	      while(abs(En.at(i)) > min && abs(En.at(i)) < max && (i+1) < lmax)
		{
		  ++ni;
		  ++i;
		  min = abs(En.at(i+1)) - 0.1;
		  max = abs(En.at(i+1)) + 0.1;
		}
	    }
	  ++i;
	  //alpha = alpha * (Kd - ni)/Kd;
	  alpha = alpha * Kd/(Kd+ni);
	  std::cout<<alpha<<std::endl;
	  //alpha = Kd / (Kd + n);                                             
	  list_X.push_back(alpha);
	  mi += ni;
	  //list_N.insert(list_N.end(), mi);  //to monitor shrincage ratio of the posterior  
	}
      En = Ena;
      i =0;
    }
  else
    {
      if(P > 0)
	{
	  for (i=0; i < imax; ++i)
	    {
	      alpha = alpha * (1 - P/(Kd+1));
	      list_X.insert(list_X.end(), alpha);
	    }
	}
      else
	{
	  for (i=0; i < imax; ++i)
	    {
	      alpha = alpha * (Kd/(Kd+1));
	      list_X.insert(list_X.end(), alpha);
	    }
	}
    }
  
  list_X.insert(list_X.begin(), 1);
  list_X.insert(list_X.end(), 0);
  
  std::cout<<"size of list_X "<<list_X.size()<<" and En "<<En.size()<<std::endl;
  //end list_X
  
  std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> > (end - start);
  
  std::cout<<"It took "<< time_span.count() << "seconds for the whole program to run"<<std::endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  //Renormalisation constant
  
  double Li = 0;
  i = 0;
  for(itp = En.begin(); itp != En.end(); ++itp)
    {
      Li +=  0.5 * ( list_X.at(i) - list_X.at(i+2) ) ;
      ++i;
    }                                                                                                           
  
  //std::cout<<"Li: "<<Li<<std::endl;                                                                            
  //int S = list_X.size();
  double ground = Li;//0.99*100000*list_X.at(S-5);
  std::cout<<"ground"<<ground<<std::endl;
  double renorm = 1;//1/ground;
  //end of renormalisation
  
  //analysis 
  int normal = 1; //change it to whatever you like
  double begin_loop = 1/T_init;
  double end_loop = 1/T_fin;
  double step_loop = begin_loop - 1/(T_init + T_step);
  std::vector<double> density;
  std::vector<double> density_t(imax,0);
  std::vector<double> capacity;
  std::vector<double> capacity_t;
  std::vector<double>::iterator it;
  std::vector<double>::iterator ite;
  double w, ct, E;
  double wt = 0;
  double Emax;
  int l;
  
  if(L==1)
    {
      Emax = 50;
      ct = n/(4*powl(Emax,0.5*n));
      i = 0;
      //std::cout<<"ct: "<<ct<<std::endl;
      
      wt = ct * (powl(Emax,0.5*n)-powl(En.at(1),0.5*n-1));
      density_t.at(i) = wt;
      ++i;
      
      for(itp = En.begin(); itp < En.end()-2; ++itp)
	{ 
	  wt = ct * (powl(*itp,0.5*n)-powl(*(itp+2),0.5*n));
	  density_t.at(i) = wt;
	  ++i;
	}
      wt = ct * (powl(*(En.end()-2),0.5*n));
      density_t.at(i) = wt;
    }
 
  output.open(dos);
  output<<"# LJ"<<L<<" Nested-Sampling by Stefano Martiniani"<<std::endl;
  output<<"Resolution: "<<K<<std::endl;
  output<<"# time taken for the plot: "<<time_span.count()<<" seconds"<<std::endl;
  output<<"# X \t ln(Likelyhood)"<<std::endl;

  i = 0;
  for(itp = En.begin(); itp < En.end(); ++itp)
    {
      w = 0.5 * ( list_X.at(i) - list_X.at(i+2)) * renorm;
      std::cout<<"w_"<<i<<" is "<<w<<std::endl; 
      density.push_back(w);
      output<<w<<"\t"<<(*itp)/normal<<"\t"<<density_t.at(i)<<std::endl;
      ++i;
    }

  output.close();
  
  ///////////////renormalise wrt ground state//////////////
  std::cout<<"En.back()"<<En.back()<<std::endl;
  for(ite = En.begin(); ite!= En.end(); ++ite)
  {
    *ite = *ite - En.back();
  }
  /////////////////////////////////////////////////////
  
  output.open(internal_energy);
  output<<"# LJ"<<L<<" Nested-Sampling by Stefano Martiniani"<<std::endl;
  output<<"Resolution: "<<K<<std::endl;
  output<<"# time taken for the plot: "<<time_span.count()<<" seconds"<<std::endl;
  output<<"# kT \t U"<<std::endl;
  
  std::vector<double> internal_U;
  std::vector<double> partition_foo;
  std::vector<double> free_en;
  std::vector<double> internal_U_t;
  std::vector<double> partition_foo_t;
  std::vector<double> free_en_t;
  std::vector<double> canonical_dist;
  
  double sum_E, sum_Et;
  double sum_Z, sum_Zt;
  double c; //must go down to 0.01, c = 1/kT
  double gc = 0;
  double g = 0;
  if (L==1)
    {
      gc = 1;//1/( pow(M_PI,0.5*n) * std::tgamma(0.5*n));    
    }
  double U,Ut;
  
  l=0;
  for(c = begin_loop; c > end_loop; c = c-step_loop)
    {
      sum_E = 0; 
      sum_Z = 0;
      sum_Et = 0; 
      sum_Zt = 0;
      l = 0;
      for(ite = En.begin(); ite != En.end(); ++ite)
        {
          E = *ite;
	  g = density_t.at(l);//;powl(E,(0.5*n - 1));
	  sum_E += E * density.at(l) * exp(-E*c);
          sum_Z += density.at(l) * exp(-E*c);
          sum_Et += E * gc * g * exp(-E*c);
          sum_Zt += gc * g * exp(-E*c);
	  ++l;
        }
      U = sum_E / sum_Z;
      Ut = sum_Et / sum_Zt;
      internal_U.push_back(U);
      internal_U_t.push_back(Ut);
      partition_foo.push_back(sum_Z);
      partition_foo_t.push_back(sum_Zt);
      output<<(1/c)<<"\t"<<(U/normal)<<"\t"<<(Ut/normal)<<std::endl;
    }

  output.close();

  output.open(heat_capacity);
  output<<"# LJ"<<L<<" Nested-Sampling by Stefano Martiniani"<<std::endl;
  output<<"Resolution: "<<K<<std::endl;
  output<<"# time taken for the plot: "<<time_span.count()<<" seconds"<<std::endl;
  output<<"# kT \t Cv/N"<<std::endl;

  double Cv, Cvt, Tc, kin_e;
  
  if (L > 1)
    {
      kin_e = (3*L - 6) / 2;
    }
  else
    {
      kin_e = 0;
    }

  i = 0;
  for(c = begin_loop; c > end_loop; c = c-step_loop)
    {
      sum_E = 0;
      sum_Et = 0;
      sum_Z = 0; 
      sum_Zt = 0;
      l = 0;
      for(ite = En.begin(); ite != En.end(); ++ite)
        {
          E = *ite;
	  g = density_t.at(l);//powl(E,0.5*n - 1);
          sum_E += (E * E * density.at(l) * exp(-E*c));
          sum_Et += gc * E * E * g * exp(-E*c);
	  ++l;
        }
      U = internal_U.at(i);
      Ut = internal_U_t.at(i);
      sum_Z = partition_foo.at(i);
      sum_Zt = partition_foo_t.at(i);
      Cv = ((sum_E/sum_Z) - U*U)*c*c + kin_e;
      Cvt = ((sum_Et/sum_Zt) - Ut*Ut)*c*c + kin_e;
      capacity.push_back(Cv);
      capacity_t.push_back(Cvt);      
      output<<(1/c)<<"\t"<<(Cv/normal)<<"\t"<<(Cvt/normal)<<std::endl;
      ++i;
    }

  output.close();

  it = std::max_element(capacity.begin(), capacity.end());
  i = it - capacity.begin();
  Tc = 1/(begin_loop - (i*step_loop));
  std::cout<<"Tc: "<<Tc<<std::endl;

  output.open(free_energy);
  output<<"# LJ"<<L<<" Nested-Sampling by Stefano Martiniani"<<std::endl;
  output<<"Resolution: "<<K<<std::endl;
  output<<"# time taken for the plot: "<<time_span.count()<<" seconds"<<std::endl;
  output<<"# kT \t F/N"<<std::endl;

  double F, Ft;
  i = 0;
  for(c = begin_loop; c > end_loop; c = c-step_loop)
    {
      sum_Z = partition_foo.at(i);
      sum_Zt = partition_foo_t.at(i);
      F = -(1/c)*( log(sum_Z));
      Ft = -(1/c)*( log(sum_Zt));
      free_en.push_back(F);
      free_en_t.push_back(Ft);

      output<<(1/c)<<"\t"<<(F/normal)<<"\t"<<(Ft/normal)<<std::endl;
      ++i;
    }

  output.close();

  output.open(entropy);
  output<<"# LJ"<<L<<" Nested-Sampling by Stefano Martiniani"<<std::endl;
  output<<"Resolution: "<<K<<std::endl;
  output<<"# time taken for the plot: "<<time_span.count()<<" seconds"<<std::endl;
  output<<"# kT \t S/N"<<std::endl;

  double Se,Set;
  i = 0;
  for(c = begin_loop; c > end_loop; c = c-step_loop)
    {
      U = internal_U.at(i);
      F = free_en.at(i);
      Ut = internal_U_t.at(i);
      Ft = free_en_t.at(i);
      Se = (U - F)*c;
      Set = (Ut - Ft)*c;
      output<<(1/c)<<"\t"<<(Se/normal)<<"\t"<<(Set/normal)<<std::endl;
      ++i;
    }

  output.close();
  /*
  output.open(canonical_distribution2);
  output<<"# LJ"<<L<<" Nested-Sampling by Stefano Martiniani"<<std::endl;
  output<<"Resolution: "<<K<<std::endl;
  output<<"# time taken for the plot: "<<time_span.count()<<" seconds"<<std::endl;
  output<<"# E/N \t D"<<std::endl;

  double D;
  double D_renormal;

  i = 0;
  for(ite = En.begin(); ite != En.end(); ++ite)
    {
      E = *ite;
      D = density.at(i) * exp(-E/Tc);
      canonical_dist.push_back(D);
      ++i;
    }

  D_renormal = *std::max_element(canonical_dist.begin(), canonical_dist.end());

  i = 0;
  for(ite = En.begin(); ite != En.end(); ++ite)
    {
      E = *ite;
      D = canonical_dist.at(i)/D_renormal;
      output<<(E/normal)<<"\t"<<D<<std::endl;
      ++i;
    }

  output.close();

  canonical_dist.clear();

  output.open(canonical_distribution22);
  output<<"# LJ"<<L<<" Nested-Sampling by Stefano Martiniani"<<std::endl;
  output<<"Resolution: "<<K<<std::endl;
  output<<"# time taken for the plot: "<<time_span.count()<<" seconds"<<std::endl;
  output<<"# E/N \t D"<<std::endl;

  i = 0;
  for(ite = En.begin(); ite != En.end(); ++ite)
    {
      E = *ite;
      D = density.at(i) * exp(-E/(Tc*0.99));
      canonical_dist.push_back(D);
      ++i;
    }

  D_renormal = *std::max_element(canonical_dist.begin(), canonical_dist.end());

  i = 0;
  for(ite = En.begin(); ite != En.end(); ++ite)
    {
      E = *ite;
      D = canonical_dist.at(i)/D_renormal;
      output<<(E/normal)<<"\t"<<D<<std::endl;
      ++i;
    }

  output.close();

  canonical_dist.clear();
 
  output.open(canonical_distribution24);
  output<<"# LJ"<<L<<" Nested-Sampling by Stefano Martiniani"<<std::endl;
  output<<"Resolution: "<<K<<std::endl;
  output<<"# time taken for the plot: "<<time_span.count()<<" seconds"<<std::endl;
  output<<"# E/N \t D"<<std::endl;

  i = 0;
  for(ite = En.begin(); ite != En.end(); ++ite)
    {
      E = *ite;
      D = density.at(i) * exp(-E/Tc*1.01);
      canonical_dist.push_back(D);
      ++i;
    }

  D_renormal = *std::max_element(canonical_dist.begin(), canonical_dist.end());

  i = 0;
  for(ite = En.begin(); ite != En.end(); ++ite)
    {
      E = *ite;
      D = canonical_dist.at(i)/D_renormal;
      output<<(E/normal)<<"\t"<<D<<std::endl;
      ++i;
    }

  output.close();
  */
  canonical_dist.clear();
  
  En.clear();
  list_X.clear();
    
  return 0;
}
