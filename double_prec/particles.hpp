#ifndef __PARTICLES_HPP
#define __PARTICLES_HPP

/********************************************//**
 * Class that stores particle positions, velocities
 * and forces in 1D dynamically allocated arrays.
 ***********************************************/

#include <stdio.h>
#include <cmath>
using namespace std;

#define MX_OVRLP 0.65


class Particles {
private:
  int m_N; // Number of particles 
  double m_R; // Raius of the particles 
  double m_M; // Mass of the partilces 
  double m_Q; // Charge of the particles 

  double* m_pos; //Position of the particles (linearized 3*N array)
  double* m_mom; //Momentum of the particles (linearized 3*N array)
  double* m_force; //Force on the particles (linearized 3*N array)

public:
  Particles();
  Particles(int N, double R, double M, double Q);
  ~Particles();

  //Getter methods
  int getNoOfParticles();
  int getNoOfParticles() const;
  double getRadius();
  double getRadius() const;
  double getMass();
  double getMass() const;
  double getCharge();
  double getCharge() const;

  double getPosition(int i, int dim);
  double getPosition(int i, int dim) const;
  double getMomentum(int i, int dim);
  double getMomentum(int i, int dim) const;
  double getForce(int i, int dim);
  double getForce(int i, int dim) const;

  double* getPositions();
  double* getPositions() const;
  double* getMomenta();
  double* getMomenta() const;
  double* getForces();
  double* getForces() const;

  //Setter methods
  void setNoOfParticles(int N);
  void setRadius(double R);
  void setMass(double M);
  void setCharge(double Q);
  
  void setPosition(int i, int dim, double val);
  void setMomentum(int i, int dim, double val);
  void setForce(int i, int dim, double val);

  void setForces(double* F, int size);
  
  //Other methods
  void addToForce(int i, int dim, double val);
  void setForcesToZero();
};

#endif
