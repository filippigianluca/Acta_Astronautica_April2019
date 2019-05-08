#ifndef SPHERICALSHAPING
#define SPHERICALSHAPING


#include <Eigen/Core>
#include <cmath>
#include <map>
#include "loaddata.h"


//---------------------------------------------------------------------//
// Spherical Shaping Function
// Return the variable isLoopConverged
// The first inputs are passed by reference: they will be the outputs; the input with const are the real inputs

int SphericalShapingFunction(
        Eigen::ArrayXd& thetaVec, Eigen::ArrayXd& R, Eigen::ArrayXd& Phi,
        Eigen::MatrixXd& u_mat, double& DV, int& cont,
        bodies DepBody, bodies ArrBody, const double t_dep, const double TOF,
        const int nr, const double toll, const double thetaStep, const int maxIter,
        double a2, const double h);

#endif // SPHERICALSHAPING

