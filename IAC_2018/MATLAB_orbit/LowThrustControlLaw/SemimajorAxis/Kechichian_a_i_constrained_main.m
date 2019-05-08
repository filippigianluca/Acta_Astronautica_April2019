% =========================================================================
% Thrusting law for minimum time transfer with variation of semimajor axis
% and inclination with contraint on the maximum semimajor axis.
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
a_initial = 6563;
% Final semimajor axis [km[
a_final   = 6878;
% Limit semimajor axis [km]
a_limit = 7000;

% Initial inclination 
i_initial = 0 * pi/180;
% Final inclination
i_final =  69.607* pi/180;

% Low thrust acceleration [km/s^2]
f = 3.5e-8;

% Plot?
plot_flag = 1;


% =========================================================================
% Constants
% =========================================================================
% Gravitational constant
mu = astro_constants(13);


% Delta_V [km/s] and ToF [h] required for the transfer
[DeltaV_equation, ToF] = Kechichian_a_i_constrained(a_initial, a_final, a_limit, i_initial, i_final, f, mu, plot_flag);