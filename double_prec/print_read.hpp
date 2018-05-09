#ifndef __PRINT_READ_HPP
#define __PRINT_READ_HPP


#include <fstream>
#include <iomanip>
#include "particles.hpp"

/**
 *Functions used for input/output from/to files
 */

bool print_conf(const Particles& part, double L, bool fin);

bool chrg_conf(Particles& part, double& L, double q);

bool print_out(double Time, double telap, long tcpu, int iprint, double U, double K, double inU, double inK, double preU, double preK);

#endif
