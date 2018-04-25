#include "eq_of_motion.hpp"



void moveX(Particles& part, double Dt){
  int N = part.getNoOfParticles();
  double mass = part.getMass();
  
  for(int i = 0; i < N; ++i){      
    for(int a = 0; a < 3; a++){
      double x = part.getPosition(i, a) + Dt * part.getMomentum(i, a) / mass; 
      part.setPosition(i, a, x);
    }
  }
}


void moveP(Particles& part, double Dt, double& K){
  int N = part.getNoOfParticles();
  double mass = part.getMass();
  
  double Dt2 = Dt / 2;
  double mom_sq = 0.;

  for(int i = 0; i < N; ++i){ 
    double psq = 0.;
    for(int a = 0; a < 3; a++){
      double p = part.getMomentum(i, a) + Dt2 * part.getForce(i, a); 
      part.setMomentum(i, a, p);
      psq += p * p;
    }
    mom_sq += psq; 
  }
  K = mom_sq / (2*mass);
}
