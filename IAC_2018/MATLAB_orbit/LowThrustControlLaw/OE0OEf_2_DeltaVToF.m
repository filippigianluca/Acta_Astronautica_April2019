% =========================================================================
% DeltaV and ToF required for variation of orbital elements
% Inputs: initial orbital elements, final orbital elements
% Outputs: ToF, DeltaV to change each orbital element separately
% Model: low thrust laws from the literature
% =========================================================================
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

clear all


addpath('SemimajorAxis')
addpath('Eccentricity')
addpath('Inclination')
addpath('ArgumentOfPerigee')
addpath('RightAscensionAscendingNode')


%% =========================================================================
% Constants:
% =========================================================================
% Gravitational acceleration [km^3/s^2]
constants.mu_dim = 398600;

% Adimensional gravitational acceleration
constants.mu = 1;

% Time constants [s]
constants.TU = 806.78;

% Length constant [km]
constants.DU = 6378.136;

% Gravity acceleration [m/s^2]
constants.g0_m_s = 9.8;

% Gravity acceleration [km/s^2]
constants.g0_km_s = constants.g0_m_s * 1e-3;

% Gravity acceleration [DU/TU^2]
constants.g0 = constants.g0_m_s * 1e-3 * constants.TU^2 / constants.DU;

% No J2 considered
constants.J2 = 0;

% Earth radius [km]
constants.R_Earth_dim = 6378;
constants.R_Earth = 1;


%% =========================================================================
% Input
% =========================================================================

% Spacecraft mass [kg]
m0 = 1000;

% Spacecraft thrust [N]
engine.T = 0.1;

% Acceleration [m/s^2]
engine.f = engine.T / m0;

% Acceleration [km/s^2]
engine.f = engine.f * 1e-3;

% Specific impulse [s]
engine.Isp = 1600;

% Initial orbital elements [km, rad]
initial_orbit.a      = 7000;
initial_orbit.e      = 0.2;
initial_orbit.i      = 30*pi/180;
initial_orbit.omega  = 60*pi/180;
initial_orbit.Omega  = 45*pi/180;
initial_orbit.E      = 180*pi/180;

% Final orbit elements [km, rad]
final_orbit.a     = 6800;
final_orbit.e     = 0.3;
final_orbit.i     = 32*pi/180;
final_orbit.omega = 70*pi/180;
final_orbit.Omega = 50*pi/180;
final_orbit.E     = 180*pi/180;


%% Semimajor axis

[DeltaV_a, ToF_a] = a0af_2_DeltaVToF(initial_orbit,...
    final_orbit,...
    engine, constants);



%% Eccentricity

[DeltaV_ecc, ToF_ecc] = e0ef_2_DeltaVToF(initial_orbit, ...
    final_orbit,...
    engine, constants);



%% Inclination

[DeltaV_incl, ToF_incl] = i0if_2_DeltaVToF(initial_orbit,...
    final_orbit,...
    engine, constants);

%% Argument of the perigee

[DeltaV_omega, ToF_omega] = omega0omegaf_2_DeltaVToF(initial_orbit, final_orbit,...
    engine, constants);


%% Right ascension of the ascending node

[Delta_V_RAAN, ToF_RAAN] = RAAN0RAANf_2_DeltaVToF(initial_orbit, final_orbit,...
    engine, constants);

