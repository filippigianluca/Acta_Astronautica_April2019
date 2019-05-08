function [DeltaV, ToF, beta] = Pollard_e_i(initial_orbit, final_orbit, f, constants)

% Computation of DeltaV, time of flight and elevation angle for transfer
% with simultanous change of eccentricity and inclination
%
% Input: initial_orbit -> structure with initial orbital elements
%        final_orbit   -> structure with final orbital elements
%        f             -> low thrust acceleration
%        constants     ->

% Reference: Simplified Analysis of Low Thrust Orbital Manuevers, Pollard

% Marilena Di Carlo
% marilena.di-carlo@strath.ac.uk

% Assumption: continuous thrust all over the orbit (ref to Pollard)
alpha = pi/2;

% Equation (9) in Pollard for the computation of the elevation angle 
sin_beta = (final_orbit.i - initial_orbit.i) * (3 * alpha + cos(alpha) * sin(alpha));

cos_beta = 2 * cos(initial_orbit.omega) * sin(alpha) * ...
           ( log( (final_orbit.e + 1) * (initial_orbit.e -1)/((final_orbit.e - 1) * (initial_orbit.e + 1))) - final_orbit.e + initial_orbit.e);

beta = atan2(abs(sin_beta), abs(cos_beta));
beta = atan(sin_beta/cos_beta);

% Equation (5) in Pollard for the computation of the DeltaV
DeltaV = sqrt(constants.mu_dim / initial_orbit.a) * (2 * alpha * abs(asin(initial_orbit.e) - asin(final_orbit.e))) / ...
        (cos(abs(beta)) * (3 * alpha + cos(alpha) * sin(alpha)));

% Time of flight    
ToF = DeltaV / f;

% keyboard