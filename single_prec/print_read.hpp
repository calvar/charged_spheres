#ifndef __PRINT_READ_HPP
#define __PRINT_READ_HPP


#include <fstream>
#include <iomanip>
#include "particles.hpp"

/**
 *Functions used for input/output from/to files
 */

bool print_conf(const Particles& part, float L, bool fin);

bool chrg_conf(Particles& part, float& L, float q);

bool print_out(float Time, float telap, long tcpu, int iprint, float U, float K, float inU, float inK, float preU, float preK);

#endif
