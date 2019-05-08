#ifndef LOADDATA
#define LOADDATA

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
//#include "approximatePlanets.h"
//#include "ephemerisBase.h"
//#include "getEphemeris.h"

//using namespace tudat_course::sample_return_mission;
using namespace tudat::basic_astrodynamics::orbital_element_conversions;
using namespace tudat;

//Load parameters, constants, and data
// Time in day, distances in AU
const double d2r = tudat::basic_mathematics::mathematical_constants::PI/180;
const double r2d = 1.0/d2r;
const double AU2km = tudat::basic_astrodynamics::physical_constants::ASTRONOMICAL_UNIT*1e-3;
const double km2AU = 1.0/AU2km;
const double T_sid = tudat::basic_astrodynamics::physical_constants::SIDEREAL_DAY;
const double mu_S = 1.327178*1e11*pow(km2AU,3)*pow(T_sid,2);
const double t0InMJD2000 = 0.0;
const double JulDay = tudat::physical_constants::JULIAN_DAY;
const double year = 365.25;

// Define struct containing Bodies data
struct bodies{
 std::string name;
 double a;
 double e;
 double i;
 double Omega;
 double w;
 double M0;
 double theta0;
 double n;
 double t0;
};

// Select ephemeris available for Planets
//ephemeris::PlanetIndices ephemerisEarth = ephemeris::Earth3D;
//Eigen::VectorXd earthKeplerianElements = ephemeris::getKeplerianElementsOfPlanet(
//            ephemerisEarth, t0InMJD2000 );

//Earth
//body[0].name="Earth";
//body[0].a = earthKeplerianElements(semiMajorAxisIndex)*1e-3*km2AU;
//body[0].e = earthKeplerianElements(eccentricityIndex);
//body[0].i = earthKeplerianElements(inclinationIndex);
//body[0].Omega = earthKeplerianElements(longitudeOfAscendingNodeIndex);
//body[0].w = earthKeplerianElements(argumentOfPeriapsisIndex);
//body[0].theta0 = earthKeplerianElements(trueAnomalyIndex);
//body[0].M0 = convertEccentricAnomalyToMeanAnomaly( convertTrueAnomalyToEccentricAnomaly( body[0].theta0,
//                                                                                         body[0].e ), body[0].e);
//body[0].n = sqrt(mu_S/pow(body[0].a,3));


//Mars
//const ephemeris::PlanetIndices ephemerisMars = ephemeris::Mars3D;
//const Eigen::VectorXd marsKeplerianElements = ephemeris::getKeplerianElementsOfPlanet(
//            ephemerisMars, t0InMJD2000 );

//body[1].name = "Mars";
//body[1].a = marsKeplerianElements(semiMajorAxisIndex)*1e-3*km2AU;
//body[1].e = marsKeplerianElements(eccentricityIndex);
//body[1].i = marsKeplerianElements(inclinationIndex);
//body[1].Omega = marsKeplerianElements(longitudeOfAscendingNodeIndex);
//body[1].w = marsKeplerianElements(argumentOfPeriapsisIndex);
//body[1].theta0 = marsKeplerianElements(trueAnomalyIndex);
//body[1].M0 = convertEccentricAnomalyToMeanAnomaly( convertTrueAnomalyToEccentricAnomaly( body[1].theta0,
//                                                                                         body[1].e ), body[1].e);
//body[1].n = sqrt(mu_S/pow(body[1].a,3));

//1989ML
//body[2].name = "1989ML";
//body[2].a = 1.272;
//body[2].e = 0.137;
//body[2].i = 4.378*d2r;
//body[2].Omega = 104.411*d2r;
//body[2].w = 183.267*d2r;
//body[2].M0 = 108.286*d2r;
//body[2].theta0 = 122.256*d2r;
//body[2].n = 6.865*1e-1*d2r;

//Tempel-1
//body[3].name = "Tempel-1";
//body[3].a = 3.124;
//body[3].e = 0.517;
//body[3].i = 10.527*d2r;
//body[3].Omega = 68.933*d2r;
//body[3].w = 178.926*d2r;
//body[3].M0 = 359.71*d2r;
//body[3].theta0 = 358.93*d2r;
//body[3].n = 1.789*1e-1*d2r;

//Neptune
//const ephemeris::PlanetIndices ephemerisNeptune = ephemeris::Neptune3D;
//const Eigen::VectorXd neptuneKeplerianElements = ephemeris::getKeplerianElementsOfPlanet(
//            ephemerisNeptune, t0InMJD2000 );
//body[4].name = "Neptune";
//body[4].a = neptuneKeplerianElements(semiMajorAxisIndex)*1e-3*km2AU;
//body[4].e = neptuneKeplerianElements(eccentricityIndex);
//body[4].i = neptuneKeplerianElements(inclinationIndex);
//body[4].Omega = neptunehKeplerianElements(longitudeOfAscendingNodeIndex);
//body[4].w = neptuneKeplerianElements(argumentOfPeriapsisIndex);
//body[4].theta0 = neptuneKeplerianElements(  trueAnomalyIndex);
//body[4].M0 = convertEccentricAnomalyToMeanAnomaly( convertTrueAnomalyToEccentricAnomaly( body[4].theta0,
//                                                                                         body[4].e ), body[4].e);
//body[4].n = sqrt(mu_S/pow(body[4].a,3));

#endif // LOADDATA

