% =========================================================================
% Maximum variation of orbital elements that is possible to obtain in a given ToF
% Inputs: initial orbital elements, ToF
% Outputs: final orbital elements, DeltaV
% Model: low thrust laws from the literature
% =========================================================================
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

clear 


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
initial_orbit.e      = 0.1;
initial_orbit.i      = 30*pi/180;
initial_orbit.omega  = 60*pi/180;
initial_orbit.Omega  = 45*pi/180;
initial_orbit.E      = 180*pi/180;

% Available time of flight [s]
ToF =  50*86400;


%% Semimajor axis

[a_final_increase, a_final_decrease, DeltaV_a_increase, DeltaV_a_decrease] = ...
    a0ToF_2_afDeltaV(initial_orbit, ToF,...
    engine, constants);


%% Eccentricity

[Delta_ecc, Delta_V_ecc] = e0ToF_2_efDeltaV(initial_orbit, ToF,...
    engine, constants);



%% Inclination

[Delta_incl, Delta_V_incl] = i0ToF_2_ifDeltaV(initial_orbit, ToF,...
    engine, constants);


%% Argument of the perigee

[Delta_omega, Delta_V_omega] = omega0ToF_2_omegafDeltaV(initial_orbit, ToF,...
    engine, constants);


%% Right ascension of the ascending node

[Delta_RAAN, Delta_V_RAAN] = RAAN0ToF_2_RAANfDeltaV(initial_orbit, ToF,...
    engine, constants);

