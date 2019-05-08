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

// Implementation of function Add Excess Velocity
//---------------------------------------------------------------//

Eigen::Vector3d addExcessVelocity(
        const double Vinf, const double alpha, const double beta, bodies DepBody, const double Theta_dep)

{

    // Calculate the velocity in tangential, normal, out of plane frame
    double Vt = Vinf*cos(beta)*cos(alpha);
    double Vn = Vinf*cos(beta)*sin(alpha);
    double Vh = Vinf*sin(beta);

    // Calculate the flight path angle at departure
    double gamma = atan( (DepBody.e*sin(Theta_dep))/(1 + DepBody.e*cos(Theta_dep)) );

    // Calculate the velocity in radial and transverse reference frame
    double Vr = Vt*sin(gamma) + Vn*cos(gamma);
    double Vtr = Vt*cos(gamma) - Vn*sin(gamma);

    // Calculate the velocity in the perifocal frame
    double V_xi = Vr*cos(Theta_dep) - Vtr*sin(Theta_dep);
    double V_eta = Vr*sin(Theta_dep) + Vtr*cos(Theta_dep);

    Eigen::Vector3d V_pf;
    V_pf << V_xi, V_eta, Vh;

    // Transform the vector from perifocal frame to geocentric frame

    Eigen::Matrix3d R1p;
    R1p <<    cos(DepBody.Omega),      -sin(DepBody.Omega),      0,
              sin(DepBody.Omega),      cos(DepBody.Omega),       0,
              0, 0, 1;

    Eigen::Matrix3d R2p;
    R2p << 1, 0, 0,
           0, cos(DepBody.i), -sin(DepBody.i),
           0, sin(DepBody.i), cos(DepBody.i);

    Eigen::Matrix3d R3p;
    R3p <<    cos(DepBody.w),      -sin(DepBody.w),      0,
              sin(DepBody.w),      cos(DepBody.w),       0,
              0, 0, 1;


    //  Rotation matrix
    Eigen::Matrix3d T_pf2ge = R1p*R2p*R3p;

    // Obtain velocity in cartesian frame
    Eigen::Vector3d VinfCart = T_pf2ge*V_pf;

    return VinfCart;

}
