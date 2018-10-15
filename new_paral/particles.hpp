#ifndef __PARTICLES_HPP
#define __PARTICLES_HPP

#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <omp.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
using namespace std;

#define MX_OVRLP 0.65
#define TSTEP_MULT 2
#define TIME_MULT 3
#define TAB_SIZE 32768 //2 to the power 15


//Class that stores particle positions and properties in 
// 1D dynamicaly alocated arays, so that copying them
// to the device is easy.************************************************
class Particles {
private:
  int Nkinds; //How many kinds of particles are there
  int* Npart; //How many particles of each kind are there
  double* M; //Mass of each kind of particle
  double* R; //Radius of each kind of particle
  double* Q; //Charge of each kind of particle

  double* pos; //Position of each particle (linearized 3*N array)
  double* mom; //Momentum of each particle (linearized 3*N array)
  double* force; //force on each particle  (linearized 3*N array)
  
public:
  Particles();
  Particles(int, const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&);
  ~Particles();

  int get_Nkinds();
  int get_N(int);
  int* get_Npart();
  double get_mass(int);
  double* get_M();
  double get_rad(int);
  double* get_R();
  double get_charge(int);
  double* get_Q();
  int get_Nkinds() const;
  int get_N(int) const;
  int* get_Npart() const;
  double get_mass(int) const;
  double* get_M() const;
  double get_rad(int) const;
  double* get_R() const;
  double get_charge(int) const;
  double* get_Q() const;
  void set_charge(int, double); 

  double get_pos(int, int, int);
  double get_pos(int, int, int) const;
  void set_pos(int, int, int, double);
  double* get_X();
  double* get_X() const;
  double get_mom(int, int, int);
  double get_mom(int, int, int) const;
  void set_mom(int, int, int, double);
  double* get_P();
  double* get_P() const;
  double get_F(int, int, int);
  double get_F(int, int, int) const;
  double* get_F();
  double* get_F() const;
  void set_zero_F();
  void add_F(int, int, int, double);
  void set_F(double*, int);
  void get_kind(int, int&, int&);
  void get_kind(int, int&, int&) const;
};

#endif
