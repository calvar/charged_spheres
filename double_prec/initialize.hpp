#ifndef __INITIALIZE_HPP
#define __INITIALIZE_HPP

/**
 *Functions used to initialize the system
 */

#include <fstream>
#include <vector>
#include <gsl/gsl_rng.h>
#include "particles.hpp"

/**
 *Test for overlaps.
 */
bool overlap_test(const Particles& part, double L, double frac);

/**
 * Set initial positions
 */
int ini_pos(Particles& part, double L, double dis);

/**
 * Set initial momenta
 */
void ini_mom(Particles& part, double temp, const gsl_rng* ran);

/**
 *Box-Muller algorithm
 */
void gauss(double variance, double& rnd1, double& rnd2);


#endif
