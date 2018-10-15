#ifndef __FORCE_HPP
#define __FORCE_HPP

#include <complex>
#include "particles.hpp"

void tabul(double *Tab, double alpha);

bool force_R(Particles& part, const double L[3], const double *Tab, double& U_sr, double& W_sr, double& U_c, double& W_c, double Dt);

void wave_v(double *Kvec, const double L[3], double alpha, int K_max);

void force_K(Particles& part, const double L[3], int K_max, double alpha, const double *Kvec, double& U, double& W);

void dot(const double a[3], const double b[3], double& ans);

#endif
