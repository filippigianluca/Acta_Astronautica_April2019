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

#include <boost/test/unit_test.hpp>
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
#include <TudatCore/Basics/testMacros.h>
#include <TudatCore/InputOutput/matrixTextFileReader.h>

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
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Mathematics/Interpolators/linearInterpolator.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "ShapingFunctions/loaddata.h"
#include "ShapingFunctions/approximatePlanets.h"
#include "ShapingFunctions/ephemerisBase.h"
#include "ShapingFunctions/getEphemeris.h"
#include "ShapingFunctions/basefuncderiv.h"
#include "ShapingFunctions/sphericalshaping.h"
#include "ShapingFunctions/timefunc.h"
#include "ShapingFunctions/addexcessvelocity.h"


using namespace tudat::basic_astrodynamics::orbital_element_conversions;
using namespace tudat_course::sample_return_mission;
using namespace std;

int main( )
{
    std::cout << "SPHERICAL SHAPING METHOD APPLIED TO REPRODUCE A GTOC2 SOLUTION" << std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;

    //Write the asteroids combinations
    Eigen::Matrix<double,12,4> ASTEROIDS;

    ASTEROIDS.row(0) << 3258076, 2000060, 2000058, 2002959;
    ASTEROIDS.row(1) << 3250293, 2000149, 2000569, 2002483;
    ASTEROIDS.row(2) << 3170221, 2000574, 2000209, 2011542;
    ASTEROIDS.row(3) << 3170221, 2001990, 2000240, 2001754;
    ASTEROIDS.row(4) << 3017309, 2000443, 2000490, 2001345;
    ASTEROIDS.row(5) << 3250293, 2000027, 2000110, 2001038;
    ASTEROIDS.row(6) << 3288933, 2001707, 2000047, 2014569;
    ASTEROIDS.row(7) << 3329255, 2000232, 2000807, 2001754;
    ASTEROIDS.row(8) << 3170221, 2000043, 2000074, 2002483;
    ASTEROIDS.row(9) << 3250293, 2000149, 2000224, 2009661;
    ASTEROIDS.row(10) << 3343104, 2000169, 2000075, 2000659;
    ASTEROIDS.row(11) << 3250293, 2000443, 2000058, 2002959;

    Eigen::Matrix<double,12,4> depDate_vector;
    Eigen::Matrix<double,12,4> arrDate_vector;

    depDate_vector.row(0) << 59870, 60373, 62069, 62737;
    depDate_vector.row(1) << 62866, 63118, 64997, 65802;
    depDate_vector.row(2) << 57372, 57849, 59587, 60139;
    depDate_vector.row(3) << 59574, 60194, 61839, 62396;
    depDate_vector.row(4) << 61073, 61348, 63268, 64101;
    depDate_vector.row(5) << 58021, 58469, 60326, 60963;
    depDate_vector.row(6) << 62201, 62544, 64534, 65484;
    depDate_vector.row(7) << 59418, 59700, 61693, 62378;
    depDate_vector.row(8) << 57561, 58106, 59717, 61025;
    depDate_vector.row(9) << 58448, 58846, 61048, 62232;
    depDate_vector.row(10) << 58246, 59215, 61821, 62642;
    depDate_vector.row(11) << 58460, 58884, 60714, 62393;

    arrDate_vector.row(0) << 60283, 61979, 62647, 63196;
    arrDate_vector.row(1) << 63028, 64907, 65712, 66662;
    arrDate_vector.row(2) << 57747, 59485, 60034, 60851;
    arrDate_vector.row(3) << 60104, 61749, 62306, 63145;
    arrDate_vector.row(4) << 61258, 63178, 64011, 64761;
    arrDate_vector.row(5) << 58379, 60236, 60872, 61735;
    arrDate_vector.row(6) << 62454, 64444, 65394, 66144;
    arrDate_vector.row(7) << 59610, 61603, 62288, 63369;
    arrDate_vector.row(8) << 57987, 59627, 60935, 61764;
    arrDate_vector.row(9) << 58752, 60826, 61991, 63175;
    arrDate_vector.row(10) << 59125, 61731, 62552, 65257;
    arrDate_vector.row(11) << 58794, 60623, 62303, 63204;

    // Load file with asteroids data (ATTENTION: ANGLE IN DEGREES)
    int jj=0;
    int nAst = 0;
    Eigen::Matrix<double,910,9> data;
    data.setZero(910,9);

    ifstream myfile;
    myfile.open("C:/TuDat/TudatBundle_v3_windows/tudatBundle/tudatApplications/sphericalShaping/SphericalShapingGTOC2/ast_ephem.txt");
    if (myfile.is_open())
    {
        std::cout << "OK - File opened" << std::endl;
        while (!myfile.eof())
        {
            myfile >> data(nAst,jj);
            if ( jj==8  )
            {
                nAst = nAst + 1;
                jj = 0;
            }
            else   {  jj = jj+1;  }

        }
        std::cout << "OK - Reading from file complete" << std::endl;
    }
    else {std::cout << "ERROR - File not found" << std::endl;}
    myfile.close();

    // Load database with asteroids combinations from GTOC2 solutions and select the desired combination
    int comb = 10;
    Eigen::VectorXd AST;
    AST = ASTEROIDS.row(comb);
    Eigen::VectorXd depDate_vec;
    depDate_vec = depDate_vector.row(comb);
    Eigen::VectorXd arrDate_vec;
    arrDate_vec = arrDate_vector.row(comb);

    // Search for the desired asteroids in the database and insert in the structure ArrBody
    bodies ArrBody[4];

    for  (int i=0; i<=3; i++)
    {
        for (int j=0; j<data.rows(); j++)
        {
            if (data(j,0)==AST(i))
            {
                ArrBody[i].name = AST(i);
                ArrBody[i].a = data(j,1);
                ArrBody[i].e = data(j,2);
                ArrBody[i].i = data(j,3)*d2r;
                ArrBody[i].Omega = data(j,4)*d2r;
                ArrBody[i].w = data(j,5)*d2r;
                ArrBody[i].M0 = data(j,6)*d2r;
                ArrBody[i].n = sqrt(mu_S/pow(ArrBody[i].a,3));
                ArrBody[i].t0 = data(j,7);
                ArrBody[i].theta0 = 0;
            }
        }
    }

    //Print informations
    std::cout << " The combination selected is the number : " << comb + 1  << std::endl;

    //Define the initial parameters
    const double m_i = 1500.0;
    const double Ceff = 4000.0*9.81*1e-3;

    //Set tolerance,step and max iterations
    const double toll = 1e-3;
    const double thetaStep = 2*tudat::basic_mathematics::mathematical_constants::PI/80;
    const int maxIter = 400;

    //Set vector with number of revolutions
    Eigen::VectorXd Nr_vect(5);
    Nr_vect << 0,1,2,3,4;

    //Set coefficient for Newton-Rapshon
    double a2 = 0.0;
    const double h = 0.001;

    // Set output directory ( \ or /)
    const std::string outputDirectory = "C:/TuDat/TudatBundle_v3_windows/tudatBundle/tudatApplications/sphericalShaping/SphericalShapingGTOC2";

    //---------------------------------------------------------------------------------------//

    // Define Departure Body (Earth)
    bodies Earth;

    //Earth
    Earth.name = "Earth";
    Earth.a = 0.999988049532578;
    Earth.e = 1.671681163160*1e-2;
    Earth.i = 0.8854353079654*1e-3*d2r;
    Earth.Omega = 175.40647696473*d2r;
    Earth.w = 287.61577546182*d2r;
    Earth.M0 = 257.60683707535*d2r;
    Earth.theta0 = 258.4751*d2r;
    Earth.n = sqrt(mu_S/pow(Earth.a,3));
    Earth.t0 = 54000;

    //    //Earth
    //    const ephemeris::PlanetIndices ephemerisEarth = ephemeris::Earth3D;
    //    const Eigen::VectorXd earthKeplerianElements = ephemeris::getKeplerianElementsOfPlanet(
    //                ephemerisEarth, t0InMJD2000 );

    //    Earth.name = "Earth";
    //    Earth.a = earthKeplerianElements(semiMajorAxisIndex)*1e-3*km2AU;
    //    Earth.e = earthKeplerianElements(eccentricityIndex);
    //    Earth.i = earthKeplerianElements(inclinationIndex);
    //    Earth.Omega = earthKeplerianElements(longitudeOfAscendingNodeIndex);
    //    Earth.w = earthKeplerianElements(argumentOfPeriapsisIndex);
    //    Earth.theta0 = earthKeplerianElements(trueAnomalyIndex);
    //    Earth.M0 = convertEccentricAnomalyToMeanAnomaly( convertTrueAnomalyToEccentricAnomaly( Earth.theta0,
    //                                                                                            Earth.e ), Earth.e);
    //    Earth.n = sqrt(mu_S/pow(Earth.a,3));
    //    Earth.t0 = 54000;

    // Transform M0 in positive value
    if (Earth.M0 < 0.0)
        Earth.M0 = Earth.M0 + 2*tudat::basic_mathematics::mathematical_constants::PI;


    //---------------------------------------------------------------------------------------//
    // Reproduce a GTOC2 trajectory

    //Initialize variables
    bodies DepBody;
    Eigen::VectorXd Vinf(4);
    Vinf << 3.5, 0, 0, 0;
    Vinf = Vinf*km2AU*T_sid;
    const double alpha = 0*d2r;
    const double beta = 0*d2r;
    Eigen::VectorXd DVlegs(4);
    Eigen::VectorXd Nrlegs(4);
    std::map<double, Eigen::MatrixXd> thetalegs;
    std::map<double, Eigen::MatrixXd> Rlegs;
    std::map<double, Eigen::MatrixXd> Philegs;
    std::map<double, Eigen::MatrixXd> uT;
    std::map<double, Eigen::MatrixXd> uN;
    std::map<double, Eigen::MatrixXd> uH;
    //    Eigen::Matrix<double,4,Eigen::Dynamic>  thetalegs;
    //    Eigen::Matrix<double,4,Eigen::Dynamic>  Rlegs;
    //    Eigen::Matrix<double,4,Eigen::Dynamic>  Philegs;
    //    Eigen::ArrayXd thetaMin;
    //    int lengthTheta = 0;

    int isLoopConverged = 0;
    int oneLegNotConv = 0;
    int  convergenceForOneNr = 0;
    double TOF = 0;
    double Nr_min = 0;
    int cont_min = 0;
    double DV_min = 1e10;
    int cont = 0;
    double DV = 0;

    //Variables for the interpolator
//    std::vector<double> independentVariable;
//    std::vector<double> dependentVariableR;
//    std::vector<double> dependentVariablePhi;
//    independentVariable.push_back(0);
//    dependentVariableR.push_back(0);
//    dependentVariablePhi.push_back(0);
//std::cout << "Here OK B " << std::endl;
//    // Declate the interpolators variables (for the next legs)
//tudat::interpolators::CubicSplineInterpolatorDouble
//                cubicSplineInterpolationR(independentVariable,dependentVariableR);
//tudat::interpolators::CubicSplineInterpolatorDouble
//                cubicSplineInterpolationPhi(independentVariable,dependentVariablePhi);
//std::cout << "Here OK C " << std::endl;
//independentVariable.clear();
//dependentVariableR.clear();
//dependentVariablePhi.clear();
//std::cout << "Here OK D " << std::endl;

    // Variables for saving the data obtained from the function
    Eigen::ArrayXd theta_min;
    Eigen::ArrayXd R_min;
    Eigen::ArrayXd Phi_min;
    Eigen::MatrixXd u_mat_min;

    Eigen::ArrayXd * thetaVec;
    Eigen::ArrayXd * R;
    Eigen::ArrayXd * Phi;
    Eigen::MatrixXd * u_mat;

    // Time variables
    time_t TIME0, TIMEcomp;
    TIME0 = time(0);

    // Start the cycle
    for (int Nleg=0; Nleg<=3; Nleg++)
    {
        convergenceForOneNr = 0;
        DV_min = 1e10;
        TOF = arrDate_vec(Nleg) - depDate_vec(Nleg);

        // Select DepartureBody
        if (Nleg==0)
            DepBody = Earth;
        else
            DepBody = ArrBody[Nleg-1];

        // Cycle on Nr
        for (int k=0; k<Nr_vect.size(); k++)
        {
            thetaVec = new Eigen::ArrayXd;
            R = new Eigen::ArrayXd;
            Phi = new Eigen::ArrayXd;
            u_mat = new Eigen::MatrixXd;

            isLoopConverged = SphericalShapingFunction(*thetaVec, *R, *Phi, *u_mat, DV, cont, DepBody, ArrBody[Nleg],
                                                       depDate_vec(Nleg),TOF,Nr_vect(k),toll,thetaStep,maxIter,a2,h,Vinf(Nleg),alpha,beta);

            if (isLoopConverged==1)
            {
                if (DV < DV_min)
                {
                    convergenceForOneNr = 1;
                    DV_min = DV;
                    Nr_min = Nr_vect(k);
                    cont_min = cont;
                    R_min = *R;
                    Phi_min = *Phi;
                    theta_min = *thetaVec;
                    u_mat_min = *u_mat;
                }
            }
            delete thetaVec;
            delete R;
            delete Phi;
            delete u_mat;
        }

        if (convergenceForOneNr == 0)
        {
            std::cout << "-------------------------------------------------------------------" << std::endl;
            std::cout << "ERROR - The loop for the leg N : " << Nleg + 1 << " has not converged; no trajectory found with the desired parameters" << std::endl;
            std::cout << "The number of iterations done is: " << cont << "/"  << maxIter << std::endl;
            std::cout << "-------------------------------------------------------------------" << std::endl;
            oneLegNotConv = 1;

        }
        else
        {
            std::cout << "-------------------------------------------------------------------" << std::endl;
            std::cout << "CONVERGENCE - The loop for the leg N : " << Nleg + 1 << " has converged; a trajectory has been found" << std::endl;
            std::cout << " With a number of revolutions equal to : " << Nr_min << std::endl;
            std::cout << "The number of iterations done is: " << cont_min << "/"  << maxIter << std::endl;
            std::cout << "-------------------------------------------------------------------" << std::endl;
std::cout << " A Here OK leg N " << Nleg + 1 << std::endl;
            DVlegs(Nleg) = DV_min*AU2km/T_sid;
            Nrlegs(Nleg) = Nr_min;
            uT[Nleg] = u_mat_min.row(0)*AU2km/pow(T_sid,2);
            uN[Nleg] = u_mat_min.row(1)*AU2km/pow(T_sid,2);
            uH[Nleg] = u_mat_min.row(2)*AU2km/pow(T_sid,2);
            thetalegs[Nleg] = theta_min;
            Rlegs[Nleg] = R_min;
            Philegs[Nleg] = Phi_min;
std::cout << " B Here OK leg N " << Nleg + 1 << std::endl;
//            // Interpolation to yield all vectors of the same length, equal to the length of the first leg
//            //Check the length of theta_min, which will be different to the length of theta_min for the previous leg.
//            // Create a vector with the same length of the first to add to the matrix
//            //(the first leg determines the length of the next thetavectors)
//            if (Nleg == 0)
//            {
//                lengthTheta = theta_min.size();
//                thetalegs.row(Nleg) = theta_min;
//                Rlegs.row(Nleg) = R_min;
//                Philegs.row(Nleg) = Phi_min;
//std::cout << "OK for 1st leg" << std::endl;
//            }
//            else
//            {
//std::cout << "Before interpolation " << Nleg << std::endl;
//                //Interpolate
//                thetaMin = Eigen::ArrayXd::LinSpaced(lengthTheta,theta_min(0),theta_min(lengthTheta-1));
//                thetalegs.row(Nleg) = thetaMin;
//std::cout << "Here OK a " << std::endl;
//                for (int t = 0; t <= lengthTheta; t++ )
//                {
//                    independentVariable.push_back( theta_min(t));
//                    dependentVariableR.push_back( R_min(t) );
//                    dependentVariablePhi.push_back( Phi_min(t) );
//                }
//std::cout << "Here OK b " << std::endl;
//tudat::interpolators::LinearInterpolatorDouble
//                        cubicSplineInterpolationR(independentVariable,dependentVariableR, tudat::interpolators::huntingAlgorithm );
//tudat::interpolators::LinearInterpolatorDouble
//                        cubicSplineInterpolationPhi(independentVariable,dependentVariablePhi, tudat::interpolators::huntingAlgorithm);
//std::cout << "Here OK c " << std::endl;
//                for (int l = 0; l < lengthTheta; l++ )
//                {
//                    Rlegs(Nleg,l) = cubicSplineInterpolationR.interpolate( thetaMin(l) );
//                    Philegs(Nleg,l) = cubicSplineInterpolationPhi.interpolate( thetaMin(l)  );
//                }
//std::cout << "Here OK d " << std::endl;
//                independentVariable.clear();
//                dependentVariableR.clear();
//                dependentVariablePhi.clear();
//std::cout << "After interpolation " << Nleg << std::endl;
//            }

        }

    }

    TIMEcomp = time(0);

    // Print results and in case save data to files
    if (oneLegNotConv == 1)
    {
        std::cout << "-------------------------------------------------------------------" << std::endl;
        std::cout << "ERROR - At least one of the four legs has not converged" << std::endl;
        std::cout << "-------------------------------------------------------------------" << std::endl;
    }
    else
    {
        std::cout << "-------------------------------------------------------------------" << std::endl;
        std::cout << "OK - All the four legs have converged" << std::endl;
        std::cout << "Total computational time: " << difftime(TIMEcomp, TIME0)/60 << " min " << std::endl;
        std::cout << "-------------------------------------------------------------------" << std::endl;

        // Transform results from AU and day to km and s
        double DVtot;
        double m_f;
        DVtot = DVlegs(0) + DVlegs(1) + DVlegs(2) + DVlegs(3);
        m_f = m_i*exp(-DVtot/Ceff);

        //Print data
        std::cout << "The DeltaV required and the final mass are " << std::endl;
        std::cout << "Performance index, J:  kg/yr : " << m_f/((arrDate_vec(3) - depDate_vec(0))/year) << " kg/yr " << std::endl;
        std::cout << "DeltaV (km/s) : " << DVtot << " - Final mass (kg) : "  << m_f << std::endl;
        std::cout << "-------------------------------------------------------------------" << std::endl;

        // Write the results to file using writeDataMapToTextFile to export the std::map into
        // .dat files

        // Export the thetalegs matrix
        tudat::input_output::writeDataMapToTextFile( thetalegs,
                                                     "thetalegs.dat",
                                                     outputDirectory,
                                                     "",
                                                     15,
                                                     15,
                                                     "" );
        // Export the Rlegs matrix
        tudat::input_output::writeDataMapToTextFile( Rlegs,
                                                     "Rlegs.dat",
                                                     outputDirectory,
                                                     "",
                                                     15,
                                                     15,
                                                     "" );
        // Export the Philegs matrix
        tudat::input_output::writeDataMapToTextFile( Philegs,
                                                     "Philegs.dat",
                                                     outputDirectory,
                                                     "",
                                                     15,
                                                     15,
                                                     "" );

        // Export the uT vector
        tudat::input_output::writeDataMapToTextFile( uT,
                                                     "uT.dat",
                                                     outputDirectory,
                                                     "",
                                                     15,
                                                     15,
                                                     "" );
        // Export the uN vector
        tudat::input_output::writeDataMapToTextFile( uN,
                                                     "uN.dat",
                                                     outputDirectory,
                                                     "",
                                                     15,
                                                     15,
                                                     "" );
        // Export the uH vector
        tudat::input_output::writeDataMapToTextFile( uH,
                                                     "uH.dat",
                                                     outputDirectory,
                                                     "",
                                                     15,
                                                     15,
                                                     "" );

//        -----------------------------------------------------------------------------------//

//         Write the results to comma separated files. These are ready to be read in by MATLAB using
//         the csvread() function and the matrix can be plotted directly as a mesh.

//        // Set output format for matrix output.
//        Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );

//        // Set absolute path to file containing theta.
//        const std::string thetaAbsolutePath = outputDirectory + "thetalegs.csv";

//        // Export the deltaV matrix.
//        std::ofstream exportFile1( thetaAbsolutePath.c_str( ) );
//        exportFile1 << thetalegs.format( csvFormat );
//        exportFile1.close( );

//        // Set absolute path to file containing R.
//        const std::string RAbsolutePath = outputDirectory + "Rlegs.csv";

//        // Export the departure dates vector.
//        std::ofstream exportFile2( RAbsolutePath.c_str( ) );
//        exportFile2.precision( 15 );
//        exportFile2 << Rlegs.format( csvFormat );
//        exportFile2.close( );

//        // Set absolute path to file containing phi.
//        const std::string PhiAbsolutePath  = outputDirectory + "Philegs.csv";

//        // Export the transfer times vector.
//        std::ofstream exportFile3( PhiAbsolutePath.c_str( ) );
//        exportFile3.precision( 15 );
//        exportFile3 << Philegs.format( csvFormat );;
//        exportFile3.close( );

    }



}
