#ifndef __EQ_OF_MOTION_HPP
#define __EQ_OF_MOTION_HPP

#include "particles.hpp"

/**
 *Function that updates the positions
 */
void moveX(Particles& part, double Dt);

/**
 *Function that updates the momenta
 */
void moveP(Particles& part, double Dt, double& K);


#endif
