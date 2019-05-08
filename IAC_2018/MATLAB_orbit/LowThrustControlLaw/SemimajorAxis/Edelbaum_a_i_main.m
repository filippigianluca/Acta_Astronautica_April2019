% =========================================================================
% Thrusting law for minimum time transfer with variation of semimajor axis
% and inclination. 
% User has to define the input (initial and final orbital elements,
% low-thrst acceleration).
% The script plot the variation of the orbital elements and return the
% DeltaV required (workspace)
% =========================================================================
% References: 
% T. Edelbaum, "Propulsion requirements for controllable satellites"
% Kechichian "analytic representations of optimal low-thrust transfer in
% circular orbit"

% Marilena Di Carlo, 2015
% marilena.di-carlo@strath.ac.uk

clear
close all
clc

addpath(genpath('..\..\spaceart_toolbox'));

% =========================================================================
% input (user defined)
% =========================================================================

% Initial semimajor axis [km]
a_initial = 8000;
% Final semimajor axis [km[
a_final   = 7000;

% Initial inclination 
i_initial = 10* pi/180;
% Final inclination
i_final = 10 * pi/180;

% Low thrust acceleration [km/s^2]
f = 1e-7;

% Plot?
plot_flag = 1;

% =========================================================================
% Constants
% =========================================================================

% Gravitational constant of the Earth
mu = astro_constants(13);


% Delta_V [km/s] and ToF [h] required for the transfer
[Delta_V_equation, ToF, a] = Edelbaum_a_i(a_initial, a_final, i_initial, i_final, f, mu, plot_flag);
