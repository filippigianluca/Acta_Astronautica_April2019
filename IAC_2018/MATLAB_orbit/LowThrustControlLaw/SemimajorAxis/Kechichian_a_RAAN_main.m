% =========================================================================
% Thrusting law for minimum time transfer with variation of semimajor axis
% and right ascension.
% User has to define the input (initial and final orbital elements,
% low-thrst acceleration).
% The script plot the variation of the orbital elements and return the
% DeltaV required (workspace)
% =========================================================================
% Reference: 
% Kechichian "analytic representations of optimal low-thrust transfer in
% circular orbit"

% Marilena Di Carlo, 2015
% marilena.di-carlo@strath.ac.uk

clear
close all
addpath(genpath('..\..\spaceart_toolbox'));


% =========================================================================
% input (user defined)
% =========================================================================
% Initial semimajor axis [km]
a_initial = 9000;
% Final semimajor axis [km[
a_final   = 8000;
% Initial RAAN
RAAN_initial = 40 * pi/180;
% Final RAAN
RAAN_final = 40 * pi/180;

% Inclination
incl = 50 * pi/180;

% Low thrust acceleration [km/s^2]
f = 1.5e-7;


% plot?
plot_flag = 1;


% =========================================================================
% Constants
% =========================================================================
% Gravitational constant
mu = astro_constants(13);



% Delta_V [km/s] and ToF [h] required for the transfer
[Delta_V_equation, ToF] = Kechichian_a_RAAN(a_initial, a_final, RAAN_initial, RAAN_final, incl, f, mu, plot_flag);

