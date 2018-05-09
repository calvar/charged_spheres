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
  float m_R; // Raius of the particles 
  float m_M; // Mass of the partilces 
  float m_Q; // Charge of the particles 

  float* m_pos; //Position of the particles (linearized 3*N array)
  float* m_mom; //Momentum of the particles (linearized 3*N array)
  float* m_force; //Force on the particles (linearized 3*N array)

public:
  Particles();
  Particles(int N, float R, float M, float Q);
  ~Particles();

  //Getter methods
  int getNoOfParticles();
  int getNoOfParticles() const;
  float getRadius();
  float getRadius() const;
  float getMass();
  float getMass() const;
  float getCharge();
  float getCharge() const;

  float getPosition(int i, int dim);
  float getPosition(int i, int dim) const;
  float getMomentum(int i, int dim);
  float getMomentum(int i, int dim) const;
  float getForce(int i, int dim);
  float getForce(int i, int dim) const;

  float* getPositions();
  float* getPositions() const;
  float* getMomenta();
  float* getMomenta() const;
  float* getForces();
  float* getForces() const;

  //Setter methods
  void setNoOfParticles(int N);
  void setRadius(float R);
  void setMass(float M);
  void setCharge(float Q);
  
  void setPosition(int i, int dim, float val);
  void setMomentum(int i, int dim, float val);
  void setForce(int i, int dim, float val);

  void setForces(float* F, int size);
  
  //Other methods
  void addToForce(int i, int dim, float val);
  void setForcesToZero();
};

#endif
