#include <vector>
#include "force.hpp"
using namespace std;

//Read configuration
bool chrg_conf(Particles& part, double L[3]){
  std::string line;
  int NC = part.get_Nkinds();
  char name[30];
  sprintf(name, "conf.dat");
  std::ifstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }
  
  for(int n = 0; n < NC; ++n){ //charges
    double q;
    Con >> q;
    part.set_charge(n, q);
  }
  Con >> L[0] >> L[1] >> L[2];
      
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    for(int i = 0; i < NN; i++){
      for(int a = 0; a < 3; ++a){
	double p;
	Con >> p;
	part.set_pos(n, i, a, p); 
      }
      for(int a = 0; a < 3; ++a){
	double m;
	Con >> m;
	part.set_mom(n, i, a, m);
      }
    }
  }
  Con.close();

  return true;
}


int main(){
  int steps = 10;
  int Nkinds = 2; //how many kinds of particles
  vector<int> N; //No. of particles per kind
  N.push_back(50000);
  N.push_back(50000);
  vector<double> m; //masses
  m.push_back(1.);
  m.push_back(1.);
  vector<double> r; //radii
  r.push_back(0.5);
  r.push_back(0.5);
  vector<double> q; //charges
  q.push_back(1.);
  q.push_back(1.);
  double alpha = 0.06; //Ewald convergence parameter
  int Kmax = 5; //Max range in K space
  double Dt = 1.e-3; //Time step
  double dens = 0.01; // density
  //---------------------------------------------
  
  //Create particles
  Particles part(Nkinds, N, m, r, q);

  //total number of particles
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);

  double L[3];
  
  //charge configuration
  chrg_conf(part, L);

  //Short range---------------------------------------
  double U_sr, W_sr;
  
  //Real part of the force
  double U_c, W_c;

  double *Tab;
  Tab = new double[2*TAB_SIZE];
  for(int i = 0; i < 2*TAB_SIZE; ++i)
    Tab[i] = 0.;
  tabul(Tab, alpha);
  
  force_R(part, L, Tab, U_sr, W_sr, U_c, W_c, Dt);

  printf("Usr=%0.5f, Uc=%0.5f\n",U_sr,U_c);

  //Long range (K part)------------------------------
  double U_k, W_k;

  int Ksize = (Kmax+1)*(Kmax+1)*(Kmax+1);
  double *Kvec;
  Kvec = new double[Ksize];
  wave_v(Kvec, L, alpha, Kmax);

  force_K(part, L, Kmax, alpha, Kvec, U_k, W_k);
  

  delete[] Tab;
  delete[] Kvec;
  
  return 0;
}
  
