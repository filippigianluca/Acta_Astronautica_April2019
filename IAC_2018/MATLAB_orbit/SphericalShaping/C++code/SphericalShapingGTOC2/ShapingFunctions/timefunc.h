#ifndef TIMEFUNC
#define TIMEFUNC

#include <Eigen/Core>
#include <cmath>
#include <map>
#include "loaddata.h"


//---------------------------------------------------------------------//
// Time Func
// The output is isDpos
// The first inputs are passed by reference: they will be the outputs; the input with const are the real inputs

void timeFunc(Eigen::ArrayXd& TPrime1, Eigen::ArrayXd& TPrime2,
        Eigen::ArrayXd R, Eigen::ArrayXd RPrime,
        Eigen::ArrayXd D, Eigen::ArrayXd DPrime);

#endif // TIMEFUNC

