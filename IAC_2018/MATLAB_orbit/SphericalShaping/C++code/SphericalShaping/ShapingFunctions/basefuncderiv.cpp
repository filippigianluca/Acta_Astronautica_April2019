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

// Implementation of Base Function Deriv
//-------------------------------------------------------------------------//

int baseFuncDeriv(
        Eigen::ArrayXd& R, Eigen::ArrayXd& Phi, Eigen::ArrayXd& RPrime1,
        Eigen::ArrayXd& PhiPrime1, Eigen::ArrayXd& RPrime2, Eigen::ArrayXd& PhiPrime2,
        Eigen::ArrayXd& RPrime3, Eigen::ArrayXd& PhiPrime3, Eigen::ArrayXd& D, Eigen::ArrayXd& DPrime,
        const Eigen::ArrayXd theta, const Eigen::VectorXd param, const double a2)

{
    using namespace std;


    //Declare variables
    int isDpos = 1;
    double a0 = param(0);
    double a1 = param(1);
    double a3 = param(2);
    double a4 = param(3);
    double a5 = param(4);
    double a6 = param(5);
    double b0 = param(6);
    double b1 = param(7);
    double b2 = param(8);
    double b3 = param(9);

    int nElem = theta.size();

    Eigen::ArrayXd ones(nElem);
    ones.setOnes(nElem,1);

    Eigen::ArrayXd fact(nElem);

    //Calculate functions
    R = ones/(a0*ones + a1*theta + a2*pow(theta,2) + (a3*ones + a4*theta)*cos(theta) + (a5*ones + a6*theta)*sin(theta));
    Phi = (b0*ones+b1*theta)*cos(theta) + (b2*ones+b3*theta)*sin(theta);

    RPrime1 = -pow(R,2)*(a1*ones + 2.0*a2*theta - a3*sin(theta) + a4*(cos(theta) - theta*sin(theta)) + a5*cos(theta) + a6*(sin(theta) + theta*cos(theta)));
    PhiPrime1 = -b0*sin(theta) + b1*(cos(theta) - theta*sin(theta)) + b2*cos(theta) + b3*(sin(theta) + theta*cos(theta));

    RPrime2 = -pow(R,2)*(2.0*a2*ones - a3*cos(theta) + a4*(-2.0*sin(theta) - theta*cos(theta)) - a5*sin(theta) + a6*(2.0*cos(theta) - theta*sin(theta))) + 2.0*(pow(RPrime1,2))/R;
    PhiPrime2 = -b0*cos(theta) + b1*(-2.0*sin(theta) - theta*cos(theta)) - b2*sin(theta) + b3*(2.0*cos(theta) - theta*sin(theta));

    RPrime3 = 6.0*(pow(RPrime1,3))/(pow(R,2)) - 6.0*R*RPrime1*(2.0*a2*ones - a3*cos(theta) + a4*(-2.0*sin(theta) - theta*cos(theta)) - a5*sin(theta) + a6*(2.0*cos(theta) - theta*sin(theta))) - pow(R,2)*(a3*sin(theta) + a4*(-3.0*cos(theta) + theta*sin(theta)) - a5*cos(theta) + a6*(-3.0*sin(theta) - theta*cos(theta)));
    PhiPrime3 = b0*sin(theta) + b1*(-3.0*cos(theta) + theta*sin(theta)) - b2*cos(theta) + b3*(-3.0*sin(theta) - theta*cos(theta));

    fact = pow(PhiPrime1,2) + pow(cos(Phi),2);
    D = -RPrime2 + 2.0*pow(RPrime1,2)/R + RPrime1*PhiPrime1*(PhiPrime2 - sin(Phi)*cos(Phi))/fact + R*fact;
    DPrime = -RPrime3 + 4.0*(RPrime1*RPrime2)/R - 2.0*pow(RPrime1,3)/pow(R,2) +
            (RPrime2*PhiPrime1*(PhiPrime2 - sin(Phi)*cos(Phi)) + RPrime1*PhiPrime2*(PhiPrime2 - sin(Phi)*cos(Phi)) + RPrime1*PhiPrime1*(PhiPrime3 - PhiPrime1*pow(cos(Phi),2) + PhiPrime1*pow(sin(Phi),2)) )/fact
            - (RPrime1*PhiPrime1*(PhiPrime2 - sin(Phi)*cos(Phi))*(2.0*PhiPrime1*PhiPrime2 - 2.0*cos(Phi)*sin(Phi)*PhiPrime1))/pow(fact,2)
            + RPrime1*fact + R*(2.0*PhiPrime1*PhiPrime2 - 2.0*cos(Phi)*sin(Phi)*PhiPrime1);



    //Check on D
    for (int i=0; i<nElem; i++)
    {
        if (D(i) < 0.0)
        {
            isDpos = 0;
            break;
        }
    }

    return isDpos;


}
