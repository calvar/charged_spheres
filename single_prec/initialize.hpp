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
bool overlap_test(const Particles& part, float L, float frac);

/**
 * Set initial positions
 */
int ini_pos(Particles& part, float L, float dis);

/**
 * Set initial momenta
 */
void ini_mom(Particles& part, float temp, const gsl_rng* ran);

/**
 *Box-Muller algorithm
 */
void gauss(float variance, float& rnd1, float& rnd2);


#endif
