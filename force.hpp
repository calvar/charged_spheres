#ifndef __FORCE_HPP
#define __FORCE_HPP

#include "particles.hpp"

/**
 *Function to compute the anergy and virial of the system and the 
 * forces on all particles. 
 */

bool force_R(Particles& part, double L, double& Usr, double& Wsr, double& Uc, double& Wc);

#endif
  
