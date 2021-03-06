#include "print_read.hpp"

//Print configuration
bool print_conf(const Particles& part, double L, bool fin){
  int N = part.getNoOfParticles();
  char name[30];
  sprintf(name, "conf.dat");
  std::ofstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }
  Con << setiosflags(ios::fixed) << setprecision(12);
  Con << L << endl;
  
  for(int i = 0; i < N; ++i){
    for(int a = 0; a < 3; ++a)
      Con << part.getPosition(i, a) << " ";
    for(int a = 0; a < 3; ++a)
      Con << part.getMomentum(i, a) << " ";
    Con << endl; 
  }
  Con.close();

  return true;
}


//Read configuration
bool chrg_conf(Particles& part, double& L, double q){
  std::string line;
  int N = part.getNoOfParticles();
  char name[30];
  sprintf(name, "conf.dat");
  std::ifstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }

  part.setCharge(q);
  Con >> L;
      
  for(int i = 0; i < N; i++){
    for(int a = 0; a < 3; ++a){
      double r;
      Con >> r;
      part.setPosition(i, a, r); 
    }
    for(int a = 0; a < 3; ++a){
      double p;
      Con >> p;
      part.setMomentum(i, a, p);
    }
  }
  Con.close();
  return true;
}


bool print_out(double Time, double telap, long tcpu, int iprint, double U, double K, double inU, double inK, double preU, double preK){
  //time difference evaluation
  double t_elap = (telap / iprint) * 1.e-6; //time / step
  double t_cpu = (static_cast<double>(tcpu) / iprint) * 1.e-6; //cpu time / step

  char name[30];
  sprintf(name, "output.dat");
  ofstream Out(name, ios::app);
  if(! Out){
    Out.close();
    return false;
  }
  Out << Time << " " << t_elap << " " << t_cpu << " | ";
  Out << setiosflags(ios::fixed) << setprecision(12);
  Out << K << " " << U << " " << (U+K)-(inU+inK) << " "
      << (U+K)-(preU+preK) << endl;
  Out.close();
  return true;
}
