/********************************************//**
 * Code for simulating a system of N soft spheres 
 * with point charges at their centers.
 ***********************************************/

#include <stdlib.h>
#include "initialize.hpp"
#include "print_read.hpp"
#include "force.hpp"
#include "eq_of_motion.hpp"
#include "timer.hpp"

int main(int argc, char* argv[]){
  if(argc < 6){
    printf("Usage:%s <No. of steps> <No. of particles> <charge> <mass> <size of simulation box>\n",argv[0]);
    return -1;
  }
  
  int Nsteps = static_cast<int>(atoi(argv[1])); /**< Number of iterations */
  int N = static_cast<int>(atoi(argv[2])); /**< Number of particles */
  double q = static_cast<double>(atof(argv[3])); /**< Absolute value of the particle charge */
  double m = static_cast<double>(atof(argv[4])); /**< Mass of the particles */
  double L = static_cast<double>(atof(argv[5])); /**< Side of the simulation box */

  double r = 1.; //Particle radius
  
  int iprint = 100; //printing frequency

  double dis = 0.2; //initial distance between particles
  int Nside = static_cast<int>(ceil(pow(static_cast<double>(N), 1./3)));
  if(Nside*(2*r+dis) > L){
    printf("The side of the box needs to be larger than %f\n", Nside*(2*r+dis));
    return 1;
  }
  
  
  double V = L*L*L; //volume
  double dens = static_cast<double>(N)/V; //density
  
  //particles with even index are positive and uneven index negative
  if(N%2 != 0){
    printf("There must be an even number of particles to have electroneutrality.");
    return 1;
  }
  
  double Dt = 1.e-3; /**< Time step */
  
  //
  //printf("%d\n%d\n%f\n%f\n%f\n",Nsteps,N,q,m,L);
  //

  //initialize the particles--------------------------------------------
  Particles part(N, r, m, q);
  //initial positions
  int check = ini_pos(part, L, dis);
  switch(check){
  case -1:
    printf("Problem oppening iniPos.dat.\n");
    return -1;
    break;
  case 1:
    printf("Initial overlap.\n");
    return 1;
    break;
  default:
    printf("Initial positions set.\n");
    break;
  }

  //initialize random num. gen.
  gsl_rng* ran = gsl_rng_alloc(gsl_rng_mt19937); //mersenne twister rng
  long seed = 5;//time(NULL);
  gsl_rng_set(ran, seed);
  //initial momenta
  double temp = 1.;
  ini_mom(part, temp, ran);
  gsl_rng_free(ran);
  printf("Initial momenta set.\n");
  //-------------------------------------------------------------------

  //print initial configuration
  print_conf(part, L, false);
  
  //Initialize energy, pressure and forces-----------------------------
  double U=0., P=0.; //Potential energy and pressure
  double Usr=0, Wsr=0; //soft sphere contributions
  double Uc=0, Wc=0; //Coulomb contributions
  if(! force_R(part, L, Usr, Wsr, Uc, Wc)){
    return 1;
  }
  U = (Usr + Uc) / N;
  P = dens + (Wsr + Wc) / (3*V);
  
  double K = 0.; //Kinetic energy
  double mom_sq = 0.;
  for(int i = 0; i < N; ++i){
    double p[3];
    for(int a = 0; a < 3; ++a)
      p[a] = part.getMomentum(i, a);
    mom_sq += p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
  }
  K = mom_sq / (2*m*N);

  printf("Initial potential energy = %f\n",U);
  printf("Initial kinetic energy = %f\n",K);
  printf("Initial pressure = %f\n",P);
  //-------------------------------------------------------------------

  //initial energy
  double inU = U;
  double inK = K;
  double preU = U;
  double preK = K;
  
  //Initialize time couters------------------------------------------
  Timer total_t;
  double tot_t_acc=0.; long tot_t_cpu=0;
  Timer force_t;
  double f_t_acc=0.; long f_t_cpu=0;
  Timer movX_t;
  double x_t_acc=0.; long x_t_cpu=0;
  Timer movP_t;
  double p_t_acc=0.; long p_t_cpu=0;
  double t_elap; long t_cpu;
  //-----------------------------------------------------------------

  double time = 0.;
  //main cycle*********************************************************
  for(int step = 0; step < Nsteps; ++step){
    total_t.start();
    
    //LEAPFROG_VERLET------------------------------

    //Find intermediate velocities
    moveP(part, Dt, K);

    //Update positions
    movX_t.start();
    moveX(part, Dt);
    movX_t.stop(t_elap, t_cpu);
    x_t_acc += t_elap; x_t_cpu += t_cpu; 
    
    //Update forces
    force_t.start();
    if(! force_R(part, L, Usr, Wsr, Uc, Wc)){
      return 1;
    }
    force_t.stop(t_elap, t_cpu);
    f_t_acc += t_elap; f_t_cpu += t_cpu; 
    
    //Update velocities
    movP_t.start();
    moveP(part, Dt, K);
    movP_t.stop(t_elap, t_cpu);
    p_t_acc += t_elap; p_t_cpu += t_cpu; 
    
    //---------------------------------------------

    U = (Usr + Uc) / N;
    P = dens + (Wsr + Wc) / (3*V);
    K /= N;
    
    time += Dt;
    
    total_t.stop(t_elap, t_cpu);
    tot_t_acc += t_elap; tot_t_cpu += t_cpu;

    if(step%iprint == 0){
      print_out(time, t_elap, t_cpu, iprint, U, K, inU, inK, preK, preU);
      preU = U;
      preK = K;
    }

  }
  //*******************************************************************

  //
  ofstream Tfile("time.dat");
  Tfile << "Time per step: " << (tot_t_acc/(Nsteps))*1e-6 
  	<< " cpu: " << (static_cast<double>(tot_t_cpu)/(Nsteps))*1e-6 << endl;
  Tfile << "Time computing forces: " << (f_t_acc/(Nsteps))*1e-6 
  	<< " cpu: " << (static_cast<double>(f_t_cpu)/(Nsteps))*1e-6 << endl;
  Tfile << "Time updating positions: " << (x_t_acc/(Nsteps))*1e-6 
  	<< " cpu: " << (static_cast<double>(x_t_cpu)/(Nsteps))*1e-6 << endl;
  Tfile << "Time updating momenta: " << (x_t_acc/(Nsteps))*1e-6 
  	<< " cpu: " << (static_cast<double>(p_t_cpu)/(Nsteps))*1e-6 << endl;
  Tfile.close();
  //
  
  return 0;
}
