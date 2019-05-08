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
#include "basefuncderiv.h"
#include "timefunc.h"
#include "m2theta.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include "TudatCore/Mathematics/BasicMathematics/coordinateConversions.h"
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

// Implementation of function Spherical Shaping

int SphericalShapingFunction(
        Eigen::ArrayXd& thetaVec, Eigen::ArrayXd& R, Eigen::ArrayXd& Phi,
        Eigen::MatrixXd& u_mat, double& DV, int& cont,
        bodies DepBody, bodies ArrBody, const double t_dep, const double TOF,
        const int nr, const double toll, const double thetaStep, const int maxIter,
        double a2, const double h)

{
    using namespace std;
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    srand (time(NULL));


    //--------------------------------------------------------------------------------------------------//
    // Find position and velocity of Departure Body

    double M_dep = DepBody.M0 + DepBody.n*(t_dep-0.0);

    double Theta_dep = M2theta(M_dep, DepBody.e);

//    tudat::basic_astrodynamics::orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly
//            meanToEccentricAnomalyDep( DepBody.e, M_dep, useDefaultInitialGuess,
//                                    initialGuess );
//    double E_dep = meanToEccentricAnomalyDep.convert( );
//    double Theta_dep = tudat::basic_astrodynamics::orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(E_dep,DepBody.e);

    if (Theta_dep < 0.0)
    {
        Theta_dep = Theta_dep + tudat::basic_mathematics::mathematical_constants::PI;
    }

    // Transform the orbital elements into cartesian coordinates
    Eigen::VectorXd keplerianElements(6);
    keplerianElements( semiMajorAxisIndex ) = DepBody.a;
    keplerianElements( eccentricityIndex ) = DepBody.e;
    keplerianElements( inclinationIndex ) = DepBody.i;
    keplerianElements( argumentOfPeriapsisIndex ) = DepBody.w;
    keplerianElements( longitudeOfAscendingNodeIndex ) = DepBody.Omega;
    keplerianElements( trueAnomalyIndex ) = Theta_dep;
    
    Eigen::VectorXd depCartesianElements(6);
    depCartesianElements = tudat::basic_astrodynamics::orbital_element_conversions::
            convertKeplerianToCartesianElements( keplerianElements,
                                                 mu_S );
    Eigen::Vector3d dep_rC = depCartesianElements.segment(0,3);
    Eigen::Vector3d dep_vC = depCartesianElements.segment(3,3);
            
    // Transform the cartesian coordinates into spherical coordinates 
    Eigen::Vector3d dep_rS;
    Eigen::Vector3d dep_vS;
    
    dep_rS(0) = dep_rC.norm();
    dep_rS(1) = atan2(dep_rC(1),dep_rC(0));
    dep_rS(2) = asin(dep_rC(2)/dep_rS(0));
    double phi = dep_rS(2);
    double theta = dep_rS(1);

    Eigen::Matrix3d rotMat;
    rotMat << cos(phi)*cos(theta), cos(phi)*sin(theta),  sin(phi),
              -sin(theta),          cos(theta),          0,
             -sin(phi)*cos(theta), -sin(phi)*sin(theta), cos(phi);
    dep_vS = rotMat*dep_vC;

    // Find position and velocity of Arrival Body
    double M_arr = ArrBody.M0 + ArrBody.n*(t_dep+TOF-0);
    double Theta_arr = M2theta(M_arr, ArrBody.e);

//    tudat::basic_astrodynamics::orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly
//            meanToEccentricAnomalyArr( ArrBody.e, M_arr, useDefaultInitialGuess,
//                                    initialGuess );
//    double E_arr = meanToEccentricAnomalyArr.convert( );
//    double Theta_arr = tudat::basic_astrodynamics::orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(E_arr,ArrBody.e);

    if (Theta_arr < 0.0)
    {
        Theta_arr = Theta_arr + tudat::basic_mathematics::mathematical_constants::PI;
    }

    // Transform the orbital elements into cartesian coordinates

    keplerianElements( semiMajorAxisIndex ) = ArrBody.a;
    keplerianElements( eccentricityIndex ) = ArrBody.e;
    keplerianElements( inclinationIndex ) = ArrBody.i;
    keplerianElements( argumentOfPeriapsisIndex ) = ArrBody.w;
    keplerianElements( longitudeOfAscendingNodeIndex ) = ArrBody.Omega;
    keplerianElements( trueAnomalyIndex ) = Theta_arr;

    Eigen::VectorXd arrCartesianElements( 6 );
    arrCartesianElements = tudat::basic_astrodynamics::orbital_element_conversions::
            convertKeplerianToCartesianElements( keplerianElements,
                                                 mu_S );
    Eigen::Vector3d arr_rC = arrCartesianElements.segment(0,3);
    Eigen::Vector3d arr_vC = arrCartesianElements.segment(3,3);

    // Transform the cartesian coordinates into spherical coordinates
    Eigen::Vector3d arr_rS;
    Eigen::Vector3d arr_vS;

    arr_rS(0) = arr_rC.norm();
    arr_rS(1) = atan2(arr_rC(1),arr_rC(0));
    arr_rS(2) = asin(arr_rC(2)/arr_rS(0));
    phi = arr_rS(2);
    theta = arr_rS(1);

    rotMat << cos(phi)*cos(theta), cos(phi)*sin(theta),  sin(phi),
              -sin(theta),          cos(theta),          0,
             -sin(phi)*cos(theta), -sin(phi)*sin(theta), cos(phi);
    arr_vS = rotMat*arr_vC;

    //--------------------------------------------------------------------------------------------------//
    // Check conditions on angles


    if ((dep_rS(2) < -tudat::basic_mathematics::mathematical_constants::PI/2)
            || (dep_rS(2) > tudat::basic_mathematics::mathematical_constants::PI/2))
    {
        cout << "ERROR - Phi_dep is out of range" << endl;
    }

    if ((arr_rS(2) < -tudat::basic_mathematics::mathematical_constants::PI/2)
            || (dep_rS(2) > tudat::basic_mathematics::mathematical_constants::PI/2))
    {
        cout << "ERROR - Phi_arr is out of range" << endl;
    }

    // Check that theta_arrival is not smaller than theta_dep, to avoid the design of a retrograde orbit
    if (arr_rS(1) < dep_rS(1))
    {
        arr_rS(1) = arr_rS(1) + 2*tudat::basic_mathematics::mathematical_constants::PI;
    }
    //--------------------------------------------------------------------------------------------------//
    // Write boundary conditions

    //Initial
    double R_i = dep_rS(0);
    double Theta_i = dep_rS(1);
    double Phi_i = dep_rS(2);
    double TPrime_i = (R_i*cos(Phi_i))/dep_vS(1);
    double RPrime_i = dep_vS(0)*TPrime_i;
    double PhiPrime_i = dep_vS(2)*TPrime_i/R_i;

    //Final
    double R_f = arr_rS(0);
    double Theta_f = arr_rS(1) + 2*tudat::basic_mathematics::mathematical_constants::PI*nr;
    double Phi_f = arr_rS(2);
    double TPrime_f = (R_f*cos(Phi_f))/arr_vS(1);
    double RPrime_f = arr_vS(0)*TPrime_f;
    double PhiPrime_f = arr_vS(2)*TPrime_f/R_f;

    //Calculate alpha_i/f
    double alpha_i = - (RPrime_i*PhiPrime_i)/(pow(PhiPrime_i,2)+ pow(cos(Phi_i),2));
    double alpha_f = - (RPrime_f*PhiPrime_f)/(pow(PhiPrime_f,2)+ pow(cos(Phi_f),2));

    //Calculate C_i/f
    double C_i = -(mu_S*pow(TPrime_i,2))/pow(R_i,2) + 2*pow(RPrime_i,2)/R_i + R_i*(pow(PhiPrime_i,2) + pow(cos(Phi_i),2))
            - RPrime_i*PhiPrime_i*(sin(Phi_i)*cos(Phi_i))/(pow(PhiPrime_i,2)+ pow(cos(Phi_i),2));
    double C_f = -(mu_S*pow(TPrime_f,2))/pow(R_f,2) + 2*pow(RPrime_f,2)/R_f + R_f*(pow(PhiPrime_f,2) + pow(cos(Phi_f),2))
            - RPrime_f*PhiPrime_f*(sin(Phi_f)*cos(Phi_f))/(pow(PhiPrime_f,2)+ pow(cos(Phi_f),2));

    //--------------------------------------------------------------------------------------------------//
    // Set parameters and variables for the cycle

    //Set vector of theta steps, vector of parameters
    double nElements = ceil( (Theta_f - Theta_i)/ (thetaStep/2) ) + 1;
    int nElem;
    nElem = int(nElements);
    if (nElem % 2 == 0)
    {
        nElem = nElem + 1;
    }

    thetaVec = Eigen::ArrayXd::LinSpaced(nElem,Theta_i,Theta_f);
    double trueThetaStep = thetaVec(2)-thetaVec(0);

    R.setZero(nElem,1);
    Phi.setZero(nElem,1);

    //Set the vectors of coefficients x = [a0,a1,a3,..,a6,b0,..,b3]
    Eigen::VectorXd x = Eigen::VectorXd::Zero(10);

    //Create the vector B (which depends only upon the initial conditions
    Eigen::VectorXd B(10);
    B << 1.0/R_i, 1.0/R_f, Phi_i, Phi_f, -RPrime_i/pow(R_i,2), -RPrime_f/pow(R_f,2),
            PhiPrime_i, PhiPrime_f, C_i - 2.0*pow(RPrime_i,2)/R_i, C_f - 2.0*pow(RPrime_f,2)/R_f;

    Eigen::MatrixXd A(10,10);
    A <<    1.0, Theta_i, cos(Theta_i), Theta_i*cos(Theta_i), sin(Theta_i), Theta_i*sin(Theta_i), 0.0, 0.0, 0.0, 0.0,
            1.0, Theta_f, cos(Theta_f), Theta_f*cos(Theta_f), sin(Theta_f), Theta_f*sin(Theta_f), 0.0, 0.0, 0.0, 0.0,
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0, cos(Theta_i), Theta_i*cos(Theta_i), sin(Theta_i), Theta_i*sin(Theta_i),
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  cos(Theta_f), Theta_f*cos(Theta_f), sin(Theta_f), Theta_f*sin(Theta_f),
            0.0,  1.0,  -sin(Theta_i), (cos(Theta_i)-Theta_i*sin(Theta_i)), cos(Theta_i), (Theta_i*cos(Theta_i)+sin(Theta_i)), 0.0, 0.0, 0.0, 0.0,
            0.0,  1.0,  -sin(Theta_f), (cos(Theta_f)-Theta_f*sin(Theta_f)), cos(Theta_f), (Theta_f*cos(Theta_f)+sin(Theta_f)), 0.0, 0.0, 0.0, 0.0,
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  -sin(Theta_i), (cos(Theta_i)-Theta_i*sin(Theta_i)), cos(Theta_i), (Theta_i*cos(Theta_i)+sin(Theta_i)),
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  -sin(Theta_f), (cos(Theta_f)-Theta_f*sin(Theta_f)), cos(Theta_f), (Theta_f*cos(Theta_f)+sin(Theta_f)),
            0.0,  0.0,  pow(R_i,2)*cos(Theta_i), pow(R_i,2)*(2.0*sin(Theta_i)+Theta_i*cos(Theta_i)), pow(R_i,2)*sin(Theta_i), pow(R_i,2)*(Theta_i*sin(Theta_i)-2.0*cos(Theta_i)), -alpha_i*cos(Theta_i), -alpha_i*(2.0*sin(Theta_i)+Theta_i*cos(Theta_i)), -alpha_i*sin(Theta_i), -alpha_i*(Theta_i*sin(Theta_i)-2.0*cos(Theta_i)),
            0.0,  0.0,  pow(R_f,2)*cos(Theta_f), pow(R_f,2)*(2.0*sin(Theta_f)+Theta_f*cos(Theta_f)), pow(R_f,2)*sin(Theta_f), pow(R_f,2)*(Theta_f*sin(Theta_f)-2.0*cos(Theta_f)), -alpha_f*cos(Theta_f), -alpha_f*(2.0*sin(Theta_f)+Theta_f*cos(Theta_f)), -alpha_f*sin(Theta_f), -alpha_f*(Theta_f*sin(Theta_f)-2.0*cos(Theta_f));

    //Set other parameters
    double err = toll+1;
    int isDpos = 0;
    int isLoopConverged = 0; 
    double a20 = a2;
    int i = 0;
    cont = 0;

    // Declare the variables that are used inside the cycle

    Eigen::VectorXd Aa2(10);
    Eigen::VectorXd b(10);
    double TOFcomp = 0.0;
    double DeltaTOF = 0.0;
    double Ta2p = 0.0;
    double Ta2m = 0.0;
    double sum = 0.0;
    double dDeltaTOF = 0.0;


    Eigen::ArrayXd RPrime1(nElem);
    Eigen::ArrayXd RPrime2(nElem);
    Eigen::ArrayXd RPrime3(nElem);
    Eigen::ArrayXd PhiPrime1(nElem);
    Eigen::ArrayXd PhiPrime2(nElem);
    Eigen::ArrayXd PhiPrime3(nElem);
    Eigen::ArrayXd D(nElem);
    Eigen::ArrayXd DPrime(nElem);
    Eigen::ArrayXd TPrime1(nElem);
    Eigen::ArrayXd TPrime2(nElem);

    RPrime1.setOnes(nElem,1);
    RPrime2.setOnes(nElem,1);
    RPrime3.setOnes(nElem,1);
    PhiPrime1.setOnes(nElem,1);
    PhiPrime2.setOnes(nElem,1);
    PhiPrime3.setOnes(nElem,1);
    D.setOnes(nElem,1);
    DPrime.setOnes(nElem,1);
    TPrime1.setOnes(nElem,1);
    TPrime2.setOnes(nElem,1);

    Eigen::ArrayXd Rp(nElem);
    Eigen::ArrayXd Phip(nElem);
    Eigen::ArrayXd RPrime1p(nElem);
    Eigen::ArrayXd RPrime2p(nElem);
    Eigen::ArrayXd RPrime3p(nElem);
    Eigen::ArrayXd PhiPrime1p(nElem);
    Eigen::ArrayXd PhiPrime2p(nElem);
    Eigen::ArrayXd PhiPrime3p(nElem);
    Eigen::ArrayXd Dp(nElem);
    Eigen::ArrayXd DPrimep(nElem);
    Eigen::ArrayXd TPrime1p(nElem);
    Eigen::ArrayXd TPrime2p(nElem);

    Rp.setOnes(nElem,1);
    Phip.setOnes(nElem,1);
    RPrime1p.setOnes(nElem,1);
    RPrime2p.setOnes(nElem,1);
    RPrime3p.setOnes(nElem,1);
    PhiPrime1p.setOnes(nElem,1);
    PhiPrime2p.setOnes(nElem,1);
    PhiPrime3p.setOnes(nElem,1);
    Dp.setOnes(nElem,1);
    DPrimep.setOnes(nElem,1);
    TPrime1p.setOnes(nElem,1);
    TPrime2p.setOnes(nElem,1);
    int isDposp = 0;

    Eigen::ArrayXd Rm(nElem);
    Eigen::ArrayXd Phim(nElem);
    Eigen::ArrayXd RPrime1m(nElem);
    Eigen::ArrayXd RPrime2m(nElem);
    Eigen::ArrayXd RPrime3m(nElem);
    Eigen::ArrayXd PhiPrime1m(nElem);
    Eigen::ArrayXd PhiPrime2m(nElem);
    Eigen::ArrayXd PhiPrime3m(nElem);
    Eigen::ArrayXd Dm(nElem);
    Eigen::ArrayXd DPrimem(nElem);
    Eigen::ArrayXd TPrime1m(nElem);
    Eigen::ArrayXd TPrime2m(nElem);

    Rm.setOnes(nElem,1);
    Phim.setOnes(nElem,1);
    RPrime1m.setOnes(nElem,1);
    RPrime2m.setOnes(nElem,1);
    RPrime3m.setOnes(nElem,1);
    PhiPrime1m.setOnes(nElem,1);
    PhiPrime2m.setOnes(nElem,1);
    PhiPrime3m.setOnes(nElem,1);
    Dm.setOnes(nElem,1);
    DPrimem.setOnes(nElem,1);
    TPrime1m.setOnes(nElem,1);
    TPrime2m.setOnes(nElem,1);
    int isDposm = 0;

    //---------------------------------------------------------------------------------//

    //Start cycle
    while ( (err > toll) && (cont < maxIter) )
    {

        //Construct the vector Aa2
        Aa2 << a2*pow(Theta_i,2), a2*pow(Theta_f,2), 0.0, 0.0, a2*2.0*Theta_i, a2*2.0*Theta_f, 0.0, 0.0, -a2*2.0*pow(R_i,2), -a2*2.0*pow(R_f,2);

        //Solve for the parameters x
        b = B-Aa2;
        x = A.colPivHouseholderQr().solve(b);

        //Calculate the functions and their derivatives
        isDpos = baseFuncDeriv(R, Phi, RPrime1, PhiPrime1, RPrime2, PhiPrime2, RPrime3, PhiPrime3, D, DPrime, thetaVec, x, a2);

                if (isDpos==1)
        {
            //Calculate T' and T''
            timeFunc(TPrime1, TPrime2, R, RPrime1, D, DPrime);

            //Calculate the time of flight with a Simpson quadrature
            sum=0.0;
            for (i=2; i<nElem-2; i+=2)
            {
                sum = sum + 2*TPrime1(i);
            }
            for (i=1; i<nElem-1; i+=2)
            {
                sum = sum + 4*TPrime1(i);
            }
            sum = sum + TPrime1(0) + TPrime1(nElem-1);
            TOFcomp = (trueThetaStep/6.0)*sum;

            err = abs(TOFcomp - TOF)/TOF;
        }

        if ((err > toll)  && (isDpos==1))
        {
            DeltaTOF = TOFcomp - TOF;

            // Calculate derivative of DeltaTOF with a central difference scheme
            isDposp = baseFuncDeriv(Rp, Phip, RPrime1p, PhiPrime1p, RPrime2p, PhiPrime2p, RPrime3p, PhiPrime3p, Dp, DPrimep, thetaVec, x, a2+h);
            isDposm = baseFuncDeriv(Rm, Phim, RPrime1m, PhiPrime1m, RPrime2m, PhiPrime2m, RPrime3m, PhiPrime3m, Dm, DPrimem, thetaVec, x, a2-h);

            if  (isDposp == 1)
            {
                 timeFunc(TPrime1p, TPrime2p, Rp, RPrime1p, Dp, DPrimep);
                        sum=0.0;
                        for (i=2; i<nElem-2; i+=2)
                        {
                            sum = sum + 2*TPrime1p(i);
                        }
                        for (i=1; i<nElem-1; i+=2)
                        {
                            sum = sum + 4*TPrime1p(i);
                        }
                        sum = sum + TPrime1p(0) + TPrime1p(nElem-1);
                        Ta2p = (trueThetaStep/6.0)*sum;

             if (isDposm==1)
             {
                timeFunc(TPrime1m, TPrime2m, Rm, RPrime1m, Dm, DPrimem);
                         sum=0.0;
                         for (i=2; i<nElem-2; i+=2)
                         {
                             sum = sum + 2*TPrime1m(i);
                         }
                         for (i=1; i<nElem-1; i+=2)
                         {
                             sum = sum + 4*TPrime1m(i);
                         }
                         sum = sum + TPrime1m(0) + TPrime1m(nElem-1);
                         Ta2m = (trueThetaStep/6.0)*sum;

                         // Obtain derivative with finite approximation
                         dDeltaTOF = (Ta2p-Ta2m)/(2.0*h);

                         a2 = a2 + DeltaTOF/dDeltaTOF;

             }
             else
             {
                 dDeltaTOF = (Ta2p-TOFcomp)/(h);

                 a2 = a2 + DeltaTOF/dDeltaTOF;
             }

            }
            else
            {
                if (isDposm==1)
                {
                  timeFunc(TPrime1m, TPrime2m, Rm, RPrime1m, Dm, DPrimem);
                            sum=0.0;
                            for (i=2; i<nElem-2; i+=2)
                            {
                                sum = sum + 2*TPrime1m(i);
                            }
                            for (i=1; i<nElem-1; i+=2)
                            {
                                sum = sum + 4*TPrime1m(i);
                            }
                            sum = sum + TPrime1m(0) + TPrime1m(nElem-1);
                            Ta2m = (trueThetaStep/6.0)*sum;

                            dDeltaTOF = (TOFcomp - Ta2m)/(h);

                            a2 = a2 + DeltaTOF/dDeltaTOF;

                }
                else
                {
                    isDpos = 0;
                    a2 = (a20-h) + 2*h*(rand() % 100 + 1)/100;
                }
            }
            
        }
        else if ((err > toll) && (isDpos == 0))
        {
          a2 = (a20-h) + 2*h*(rand() % 100 + 1)/100;
        }
        else if ((err <= toll) && (isDpos == 1))
        {
          isLoopConverged = 1;
        }
        cont = cont + 1;

    }

    if (isLoopConverged == 1)
    {
        // Calculate thetaDot1 and thetaDot2 (first and second derivative)
        Eigen::ArrayXd thetaDot1(nElem);
        Eigen::ArrayXd thetaDot2(nElem);
        Eigen::ArrayXd ones(nElem);
        ones.setOnes(nElem,1);
        Eigen::ArrayXd zeros(nElem);
        zeros.setZero(nElem,1);

        thetaDot1 = ones/TPrime1;
        thetaDot2 = -TPrime2/(pow(TPrime1,3));

        //Compute the r,v,a vectors in the reference frame radial,orthoradial,out of plane
        Eigen::ArrayXd U(nElem);
        U = pow(PhiPrime1,2) + pow(cos(Phi),2);

        // Doesn't work with map<double, Array<double,1,Dynamic> > Vtilde_mat because then I can't take the column
        Eigen::Matrix<double,3,Eigen::Dynamic> Vtilde_mat(3,nElem);
        Vtilde_mat.row(0) = RPrime1;
        Vtilde_mat.row(1) = R*sqrt(U);
        Vtilde_mat.row(2) = zeros;

        Eigen::Matrix<double,3,Eigen::Dynamic> Atilde_mat(3,nElem);
        Atilde_mat.row(0) = RPrime2-R*U;
        Atilde_mat.row(1) = 2*RPrime1*sqrt(U)+R*PhiPrime1*(PhiPrime2-sin(Phi)*cos(Phi))/(sqrt(U));
        Atilde_mat.row(2) = (cos(Phi)*(PhiPrime2-sin(Phi)*cos(Phi)) + 2*sin(Phi)*U)*R/(sqrt(U));

        // Define the unit vectors and calculate the unit vectors in the t-n-h frame
        Eigen::Vector3d Er; Er << 1.0,0.0,0.0;
        Eigen::Vector3d Eo; Eo << 0.0,1.0,0.0;
        Eigen::Vector3d Eh; Eh << 0.0,0.0,1.0;

        Eigen::Matrix<double,3,Eigen::Dynamic> Et_mat(3,nElem);
        Eigen::Matrix<double,3,Eigen::Dynamic> En_mat(3,nElem);

        for (i=0; i<nElem; i++)
        {
            Et_mat.col(i) = Vtilde_mat.col(i)/(Vtilde_mat.col(i)).norm();
            En_mat.col(i) = Eh.cross(Et_mat.col(i));
        }

        // Compute the control acceleration (in the t-n-h frame) (in AU/day^2)
        //Eigen::Matrix<double,3,Eigen::Dynamic> u_mat(3,nElem);
        Eigen::ArrayXd uMod_vec(nElem);
        Eigen::Vector3d auxV;
        Eigen::Vector3d auxA;

        u_mat.setOnes(3,nElem);
        uMod_vec.setOnes(nElem,1);

//std::cout << "R : " << R  << std::endl;
//std::cout << "RPrime1 : " << RPrime1  << std::endl;
//std::cout << "RPrime2 : " << RPrime2  << std::endl;
//std::cout << "RPrime3 : " << RPrime3  << std::endl;
//std::cout << "Phi : " << Phi  << std::endl;
//std::cout << "PhiPrime1 : " << PhiPrime1  << std::endl;
//std::cout << "PhiPrime2 : " << PhiPrime2  << std::endl;
//std::cout << "PhiPrime3 : " << PhiPrime3  << std::endl;
//std::cout << "thetaVec " << thetaVec << std::endl;
//std::cout << "D : " << D  << std::endl;
//std::cout << "DPrime " << DPrime << std::endl;

        for (i=0; i<nElem; i++)
        {
            auxV = Vtilde_mat.col(i);
            auxA = Atilde_mat.col(i);

            u_mat(0,i) = (mu_S/pow(R(i),2))*Er.dot(Et_mat.col(i)) + thetaDot2(i)*auxV.dot(Et_mat.col(i))
                           + pow(thetaDot1(i),2)*auxA.dot(Et_mat.col(i));

            u_mat(1,i) = (mu_S/pow(R(i),2))*Er.dot(En_mat.col(i)) + pow(thetaDot1(i),2)*auxA.dot(En_mat.col(i));

            u_mat(2,i) = pow(thetaDot1(i),2)*Eh.dot(auxA);


          uMod_vec(i) = sqrt(pow(u_mat(0,i),2) + pow(u_mat(1,i),2) + pow(u_mat(2,i),2));

        }

         //Compute the DeltaV  (in AU/day)
        sum=0.0;
        for (i=2; i<nElem-2; i+=2)
        {
            sum = sum + 2*TPrime1(i)*uMod_vec(i);
        }
        for (i=1; i<nElem-1; i+=2)
        {
            sum = sum + 4*TPrime1(i)*uMod_vec(i);
        }
        sum = sum + TPrime1(0)*uMod_vec(0) + TPrime1(nElem-1)*uMod_vec(nElem-1);
        DV = (trueThetaStep/6.0)*sum;

    }
    else
    {
        R.setZero(nElem,1);
        Phi.setZero(nElem,1);
        u_mat.setZero(3,nElem);
        DV=0.0;
    }

//std::cout << " Inside func"  << std::endl;
//std::cout << "isLoopConverged : " << isLoopConverged  << std::endl;
//std::cout << "DV: " << DV  << std::endl;
    return isLoopConverged;

}
