/*SphericalShapingFunction


*/

#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <iostream>
#include <cmath>
#include <ctime>
#include <map>
#include <vector>
#include <cstdlib>

#include <boost/assign/list_of.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>
#include <Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h>
#include <Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h>
#include <Tudat/Astrodynamics/Gravitation/centralGravityModel.h>
#include <Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h>
#include <Tudat/Astrodynamics/Gravitation/thirdBodyPerturbation.h>
#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"
#include "Tudat/Astrodynamics/MissionSegments/escapeAndCapture.h"
#include <Tudat/Astrodynamics/StateDerivativeModels/cartesianStateDerivativeModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/compositeStateDerivativeModel.h>
//#include <Tudat/External/SpiceInterface/spiceEphemeris.h>
#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include <Tudat/Mathematics/Interpolators/lagrangeInterpolator.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include "ShapingFunctions/loaddata.h"
#include "ShapingFunctions/approximatePlanets.h"
#include "ShapingFunctions/ephemerisBase.h"
#include "ShapingFunctions/getEphemeris.h"
#include "ShapingFunctions/basefuncderiv.h"
#include "ShapingFunctions/sphericalshaping.h"
#include "ShapingFunctions/timefunc.h"


using namespace tudat::basic_astrodynamics::orbital_element_conversions;
using namespace tudat_course::sample_return_mission;
using namespace std;

int main( )
{

    //Define the initial parameters
    const double m_i = 1000.0;
    const double Ceff = 3000.0*9.81*1e-3;

    //Set tolerance,step and max iterations
    const double toll = 1e-3;
    const double thetaStep = 2*tudat::basic_mathematics::mathematical_constants::PI/60;
    const int maxIter = 60;

    //Set the departure time, TOF (MJD2000 and days), and num revs
    const double t_dep = 9700.0;
    const double TOF = 25500.0;
    const int nr = 1;

    //Set coefficient for Newton-Rapshon
    double a2 = 0.0;
    const double h = 0.001;


    // Set output directory ( \ or /)
    const std::string outputDirectory = "C:/TuDat/TudatBundle_v3_windows/tudatBundle/tudatApplications/sphericalShaping/SphericalShaping";

    //---------------------------------------------------------------------------------------//

    // Define bodies
    bodies body[5];

    const ephemeris::PlanetIndices ephemerisEarth = ephemeris::Earth3D;
    const Eigen::VectorXd earthKeplerianElements = ephemeris::getKeplerianElementsOfPlanet(
                ephemerisEarth, t0InMJD2000 );
    const ephemeris::PlanetIndices ephemerisMars = ephemeris::Mars3D;
    const Eigen::VectorXd marsKeplerianElements = ephemeris::getKeplerianElementsOfPlanet(
                ephemerisMars, t0InMJD2000 );
    const ephemeris::PlanetIndices ephemerisNeptune = ephemeris::Neptune3D;
    const Eigen::VectorXd neptuneKeplerianElements = ephemeris::getKeplerianElementsOfPlanet(
                ephemerisNeptune, t0InMJD2000 );


    //Earth
    body[0].name = "Earth";
    body[0].a = 1.00000011;
    body[0].e = 0.01671022;
    body[0].i = 0.000054*d2r;
    body[0].Omega = -0.0892318;
    body[0].w = 1.8857;
    body[0].M0 = 358.189*d2r;
    body[0].theta0 = 358.126*d2r;
    body[0].n = 9.85*1e-1*d2r;

    //Mars
    body[1].name = "Mars";
    body[1].a = 1.524;
    body[1].e = 0.093;
    body[1].i = 1.85*d2r;
    body[1].Omega = 49.557*d2r;
    body[1].w = 286.502*d2r;
    body[1].M0 = 19.095*d2r;
    body[1].theta0 = 23.02*d2r;
    body[1].n = 5.24*1e-1*d2r;

    //Neptune
    body[4].name = "Neptune";
    body[4].a = 30.104;
    body[4].e = 0.011;
    body[4].i = 1.768*d2r;
    body[4].Omega = 131.794*d2r;
    body[4].w = 265.647*d2r;
    body[4].M0 = 266.599*d2r;
    body[4].theta0 = 265.325*d2r;
    body[4].n = 5.969*1e-3*d2r;

//    //Earth
//    body[0].name = "Earth";
//    body[0].a = earthKeplerianElements(semiMajorAxisIndex);
//    body[0].e = earthKeplerianElements(eccentricityIndex);
//    body[0].i = earthKeplerianElements(inclinationIndex);
//    body[0].Omega = earthKeplerianElements(longitudeOfAscendingNodeIndex);
//    body[0].w = earthKeplerianElements(argumentOfPeriapsisIndex);
//    body[0].theta0 = earthKeplerianElements(trueAnomalyIndex);
//    body[0].M0 = convertEccentricAnomalyToMeanAnomaly( convertTrueAnomalyToEccentricAnomaly( body[0].theta0,
//                                                                                             body[0].e ), body[0].e);
//    body[0].n = sqrt(mu_S/pow(body[0].a,3));

//    //Mars
//    body[1].name = "Mars";
//    body[1].a = marsKeplerianElements(semiMajorAxisIndex);
//    body[1].e = marsKeplerianElements(eccentricityIndex);
//    body[1].i = marsKeplerianElements(inclinationIndex);
//    body[1].Omega = marsKeplerianElements(longitudeOfAscendingNodeIndex);
//    body[1].w = marsKeplerianElements(argumentOfPeriapsisIndex);
//    body[1].theta0 = marsKeplerianElements(trueAnomalyIndex);
//    body[1].M0 = convertEccentricAnomalyToMeanAnomaly( convertTrueAnomalyToEccentricAnomaly( body[1].theta0,
//                                                                                             body[1].e ), body[1].e);
//    body[1].n = sqrt(mu_S/pow(body[1].a,3));

//    //Neptune
//    body[4].name = "Neptune";
//    body[4].a = neptuneKeplerianElements(semiMajorAxisIndex);
//    body[4].e = neptuneKeplerianElements(eccentricityIndex);
//    body[4].i = neptuneKeplerianElements(inclinationIndex);
//    body[4].Omega = neptuneKeplerianElements(longitudeOfAscendingNodeIndex);
//    body[4].w = neptuneKeplerianElements(argumentOfPeriapsisIndex);
//    body[4].theta0 = neptuneKeplerianElements(trueAnomalyIndex);
//    body[4].M0 = convertEccentricAnomalyToMeanAnomaly( convertTrueAnomalyToEccentricAnomaly( body[4].theta0,
//                                                                                             body[4].e ), body[4].e);
//    body[4].n = sqrt(mu_S/pow(body[4].a,3));

    //1989ML
    body[2].name = "1989ML";
    body[2].a = 1.272;
    body[2].e = 0.137;
    body[2].i = 4.378*d2r;
    body[2].Omega = 104.411*d2r;
    body[2].w = 183.267*d2r;
    body[2].M0 = 108.286*d2r;
    body[2].theta0 = 122.256*d2r;
    body[2].n = 6.865*1e-1*d2r;

    //Tempel-1
    body[3].name = "Tempel-1";
    body[3].a = 3.124;
    body[3].e = 0.517;
    body[3].i = 10.527*d2r;
    body[3].Omega = 68.933*d2r;
    body[3].w = 178.926*d2r;
    body[3].M0 = 359.71*d2r;
    body[3].theta0 = 358.93*d2r;
    body[3].n = 1.789*1e-1*d2r;



    // Transform M0 in positive value
    for (int ii=0; ii<5; ii++)
    {
        if (body[ii].M0 < 0.0)
                body[ii].M0 = body[ii].M0 + 2*tudat::basic_mathematics::mathematical_constants::PI;
    }

    // Select bodies
    bodies DepBody = body[0];
    bodies ArrBody = body[4];

    //---------------------------------------------------------------------------------------//


    //Print informations
    std::cout << "SPHERICAL SHAPING METHOD TEST CASES" << std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;
    std::cout << "TRIAL ON A SINGLE TRANSFER : " << DepBody.name << "-" << ArrBody.name << std::endl;
    std::cout << "Dep Date : " << t_dep << " - TOF : " << TOF << " - N revs: " << nr << " - a2_0 : " << a2 << std::endl;

    //---------------------------------------------------------------------------------------//
    // Call the function

    //Declare variables

    int isLoopConverged = 0;
    Eigen::ArrayXd thetaVec;
    Eigen::ArrayXd R;
    Eigen::ArrayXd Phi;
    Eigen::MatrixXd u_mat;
    double DV = 0.0;
    int cont = 0;


    isLoopConverged = SphericalShapingFunction(thetaVec, R, Phi, u_mat, DV, cont, DepBody, ArrBody,
                                               t_dep,TOF,nr,toll,thetaStep,maxIter,a2,h);

    if (isLoopConverged == 0)
    {
        std::cout << "-------------------------------------------------------------------" << std::endl;
        std::cout << "ERROR - The loop has not converged; no trajectory found with the desired parameters" << std::endl;
        std::cout << "The number of iterations done is  " << cont << "/" << maxIter << std::endl;
        std::cout << "-------------------------------------------------------------------" << std::endl;
    }
    else
    {
        std::cout << "-------------------------------------------------------------------" << std::endl;
        std::cout << "CONVERGENCE - The loop has converged; a trajectory has been found" << std::endl;
        std::cout << "The number of iterations done is  " << cont << "/" << maxIter << std::endl;
        std::cout << "-------------------------------------------------------------------" << std::endl;

        // Transform results from AU and day to km and s
        DV = DV*AU2km/T_sid;
        double m_f = m_i*exp(-DV/Ceff);
        u_mat = u_mat*AU2km/pow(T_sid,2);


        //Print data
        std::cout << "The DeltaV required and the final mass are " << std::endl;
        std::cout << "DeltaV (km/s) : " << DV << " - Final mass (kg) : "  << m_f << std::endl;
        std::cout << "-------------------------------------------------------------------" << std::endl;

        // Write the results to comma separated files. These are ready to be read in by MATLAB using
        // the csvread() function and the matrix can be plotted directly as a mesh.

        // Set output format for matrix output.
        Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );

        // Set absolute path to file containing u_mat.
        const std::string u_matOutputFileAbsolutePath = outputDirectory + "u_mat.csv";

        // Export the deltaV matrix.
        std::ofstream exportFile1( u_matOutputFileAbsolutePath.c_str( ) );
        exportFile1 << u_mat.format( csvFormat );
        exportFile1.close( );


        // Set absolute path to file containing departure dates.
        const std::string thetaVecPath
                = outputDirectory + "thetaVec.csv";

        // Export the departure dates vector.
        std::ofstream exportFile4( thetaVecPath.c_str( ) );
        exportFile4.precision( 15 );
        exportFile4 << thetaVec;
        exportFile4.close( );

        // Set absolute path to file containing departure dates.
        const std::string RPath
                = outputDirectory + "R.csv";

        // Export the departure dates vector.
        std::ofstream exportFile5( RPath.c_str( ) );
        exportFile5.precision( 15 );
        exportFile5 << R;
        exportFile5.close( );

        // Set absolute path to file containing departure dates.
        const std::string PhiPath
                = outputDirectory + "Phi.csv";

        // Export the departure dates vector.
        std::ofstream exportFile6( PhiPath.c_str( ) );
        exportFile6.precision( 15 );
        exportFile6 << Phi;
        exportFile6.close( );

    }

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    // Perform grid search
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    std::cout << "-------------------------------------------------------------------" << std::endl;
    std::cout << "GRID SEARCH OF THE TRANSFER : " << DepBody.name << " - " <<  ArrBody.name <<std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;

    //Set vectors of DepTime and TOF
    Eigen::VectorXd t_dep_vect = Eigen::VectorXd::LinSpaced(51,7500,10000);
    Eigen::VectorXd TOF_vect = Eigen::VectorXd::LinSpaced(96,11000,30000);
    int Ndep = t_dep_vect.size();
    int Ntof = TOF_vect.size();

    int convPerc = 0;
    Eigen::MatrixXd DV_matrix(Ntof,Ndep);
    DV_matrix.setZero(Ntof,Ndep);

    //Parameters for the case of check on different number of revolutions
    Eigen::VectorXd Nr_vect(2);                //Have to declare a size otherwise error
    Nr_vect << 0,1;
    int Nr = Nr_vect.size();

    int convergenceForOneNr = 0;
    double DV_min = 0.0;


    //Initialise vectors that will be output
    Eigen::ArrayXd * thetaVecg;
    Eigen::ArrayXd * Rg;
    Eigen::ArrayXd * Phig;
    Eigen::MatrixXd * u_matg;
//      std::map<double, Eigen::ArrayXd> thetaVecg;
//      std::map<double, Eigen::ArrayXd> Rg;
//      std::map<double, Eigen::ArrayXd> Phig;
//      std::map<double, Eigen::MatrixXd> u_matg;
//        int comb = 0;

    // Set variables for keeping track of time
    time_t TIME0, TIMEcomp;
    TIME0 = time(0);

    // Start cycle
    for (int i=0; i<Ntof; i++)
    {
        for (int j=0; j<Ndep; j++)
        {
           DV_min = 1000000.0;
           convergenceForOneNr = 0;
           for (int k=0; k<Nr; k++)
           {
               thetaVecg = new Eigen::ArrayXd;
               Rg = new Eigen::ArrayXd;
               Phig = new Eigen::ArrayXd;
               u_matg = new Eigen::MatrixXd;

//             comb = comb + 1;
               isLoopConverged = SphericalShapingFunction(*thetaVecg, *Rg, *Phig, *u_matg, DV, cont, DepBody, ArrBody,
                                                          t_dep_vect(j),TOF_vect(i),Nr_vect(k),toll,thetaStep,maxIter,a2,h);

               delete thetaVecg;
               delete Rg;
               delete Phig;
               delete u_matg;

               if (isLoopConverged == 1)
                    {
                     if (DV < DV_min)
                        {
                         convergenceForOneNr = 1;
                         DV_min = DV;
                        }
                    }
           }
           if (convergenceForOneNr == 1)
           {
               convPerc = convPerc+1;
               DV_matrix(i,j) = DV_min;
           }
           else
           {
               DV_matrix(i,j) = TUDAT_NAN;
           }

        }
    }

    TIMEcomp = time(0);

    std::cout << "Results: " << std::endl;
    std::cout << "Total computational time: " << difftime(TIMEcomp, TIME0)/60 << " min " << std::endl;
    std::cout << "Comp time per trajectory: " << difftime(TIMEcomp, TIME0)/(Ntof*Ndep) << " s " << std::endl;
    std::cout << "Number of obtained trajectories " << convPerc << " / " << Ntof*Ndep << std::endl;
    double convD = 0;
    convD = double(convPerc);
    std::cout << "Percentage of feasible trajectories : " << convD/(Ntof*Ndep)*100 << " % " << std::endl;
    std::cout << "DeltaV of best trajectory (km/s):  " <<  DV_matrix.minCoeff()*AU2km/T_sid << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

    // Set output format for matrix output.
    Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );

    // Set absolute path to file containing deltaVs.
    const std::string deltaVOutputFileAbsolutePath = outputDirectory + "DeltaV.csv";

    // Export the deltaV matrix.
    std::ofstream exportFile7( deltaVOutputFileAbsolutePath.c_str( ) );
    exportFile7 << DV_matrix.format( csvFormat );
    exportFile7.close( );

    // Set absolute path to file containing departure dates.
    const std::string departureDateOutputFileAbsolutePath
            = outputDirectory + "DepartureDates.csv";

    // Export the departure dates vector.
    std::ofstream exportFile2( departureDateOutputFileAbsolutePath.c_str( ) );
    exportFile2.precision( 15 );
    exportFile2 << t_dep_vect.format( csvFormat );
    exportFile2.close( );

    // Set absolute path to file containing transfer times.
    const std::string transferTimesOutputFileAbsolutePath
            = outputDirectory + "TransferTimes.csv";

    // Export the transfer times vector.
    std::ofstream exportFile3( transferTimesOutputFileAbsolutePath.c_str( ) );
    exportFile3.precision( 15 );
    exportFile3 << TOF_vect.format( csvFormat );;
    exportFile3.close( );


}
