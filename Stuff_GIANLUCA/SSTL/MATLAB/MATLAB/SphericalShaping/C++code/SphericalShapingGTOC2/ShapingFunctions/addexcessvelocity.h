#ifndef ADDEXCESSVELOCITY
#define ADDEXCESSVELOCITY

#include <Eigen/Core>
#include <cmath>
#include <map>
#include "loaddata.h"


//---------------------------------------------------------------------//
// Add Excess Velocity
// The output is VinfCart
// The first inputs are passed by reference: they will be the outputs; the input with const are the real inputs
Eigen::Vector3d addExcessVelocity(
        const double Vinf, const double alpha, const double beta, bodies DepBody, const double Theta_dep);


#endif // ADDEXCESSVELOCITY

