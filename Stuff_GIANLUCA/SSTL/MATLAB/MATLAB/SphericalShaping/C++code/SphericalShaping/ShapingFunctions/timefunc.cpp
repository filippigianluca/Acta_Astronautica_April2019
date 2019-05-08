#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <cmath>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <iostream>

#include "loaddata.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include "TudatCore/Mathematics/BasicMathematics/coordinateConversions.h"
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

// Implementation of Time Func
//-------------------------------------------------------------------------//

void timeFunc(Eigen::ArrayXd& TPrime1, Eigen::ArrayXd& TPrime2,
        Eigen::ArrayXd R, Eigen::ArrayXd RPrime,
        Eigen::ArrayXd D, Eigen::ArrayXd DPrime)

{
    using namespace std;

    TPrime1 = sqrt(D*pow(R,2)/mu_S);
    TPrime2 = (2.0*D*RPrime + R*DPrime)/(2.0*sqrt(mu_S*D));

}
