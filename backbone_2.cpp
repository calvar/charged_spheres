#include <cmath>
#include <complex>
#include <random>
#include <stdio.h>
using namespace std;


///Parte REAL--------------------------------
void force_R(int N, const double* positions, double* rForce, double& sum){
  /*
    Compute the sum of all interparticle distances in 1D
  */
  sum = 0;
  for(int i = 0; i < N; ++i){
    
    //Esta forma del ciclo interno equilibra las cargas-------------------
    int mx = static_cast<int>(ceil(static_cast<float>(N-1)/2));
    if(fmod(static_cast<float>(N),2) == 0. && i >= N/2)
      mx = static_cast<int>(floor(static_cast<float>(N-1)/2));
    
    int j = i+1 - N*static_cast<int>(floor((i+1)/N + 0.5));
    int cnt = 0;
    while(cnt < mx){
      //printf("(%d,%d)\n",i,j);
      
      //Mirar que operaciones se pueden vectorizar dentro de este ciclo
      double XIJ = positions[i] - positions[j];
      sum += 1;

      rForce[i] += XIJ;
      rForce[j] -= XIJ;

      j += 1 - N*static_cast<int>(floor((j+1)/N + 0.5));
      cnt++;
    }
    //---------------------------------------------------------------------
    
  }
}



///Parte reciproca--------------------------
void force_K(int N, double Kmax, const double* positions, double* kForce, double& sum){
  /*
    compute the sum of a function of the positions
  */
  sum = 0;
  double P2 = 2.*M_PI;
  
  for(int kn = 0; kn < Kmax; ++kn){
    double Kx = P2 * kn;
    
    double K[2] = {Kx, -Kx};
    complex<double> cc(0., 0.);
    complex<double> sum = cc;

    complex<double> M[N*2];
    for(int a = 0; a < N*2; ++a){
      M[a] = cc;
    }

    //-------------------------------
    for(int i = 0; i < N; ++i){  //N es del orden de 10000+
      double xi = positions[i];
      double RIK[2];
      complex<double> z[2];
      for(int b = 0; b < 2; ++b){
	RIK[b] = xi*K[b]; 
	z[b] = polar(1., RIK[b]);
      }

      for(int b = 0; b < 2; ++b){
	M[i*2+b] += conj(z[b]);
      }
    }
    //-------------------------------
    
    //varias operaciones entre elementos de un arreglo

    //-------------------------------
    for(int i = 0; i < N; ++i){
      for(int b = 0; b < 2; ++b){
	kForce[i] += K[b] * norm(M[i*2+b]);
      }
      printf("%f\n", positions[i]);
      sum += 1;//positions[i];
    }
    //-------------------------------
    
  }
}


int main(){
  int steps = 100;
  int N = 100;
  int Kmax = 400;

  default_random_engine rng;
  uniform_real_distribution<double> distribution(0.0,1.0);

  double positions[N];
  double rForce[N];
  double kForce[N];
  for(int i = 0; i < N; ++i){
    positions[i] = 2*distribution(rng)-1;
    rForce[i] = 0.;
    rForce[i] = 0.;
  }
  
  double sumR = 0;
  force_R(N, positions, rForce, sumR);

  double sumK = 0;
  force_K(N, Kmax, positions, kForce, sumK);

  printf("%f, %f\n", sumR, sumK);
  
  return 0;
}
