#ifndef BASEFUNCDERIV
#define BASEFUNCDERIV

#include <Eigen/Core>
#include <cmath>
#include <map>
#include "loaddata.h"


//---------------------------------------------------------------------//
// Base Function Deriv
// The output is isDpos
// The first inputs are passed by reference: they will be the outputs; the input with const are the real inputs
int baseFuncDeriv(
        Eigen::ArrayXd& R, Eigen::ArrayXd& Phi, Eigen::ArrayXd& RPrime1,
        Eigen::ArrayXd& PhiPrime1, Eigen::ArrayXd& RPrime2, Eigen::ArrayXd& PhiPrime2,
        Eigen::ArrayXd& RPrime3, Eigen::ArrayXd& PhiPrime3, Eigen::ArrayXd& D, Eigen::ArrayXd& DPrime,
        const Eigen::ArrayXd theta, const Eigen::VectorXd param, const double a2);

#endif // BASEFUNCDERIV

