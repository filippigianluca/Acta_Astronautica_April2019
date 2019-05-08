function [DeltaV, ToF] = Pollard_omega(initial_orbit, final_orbit, f, constants)

% Computation of DeltaV and time of flight for transfer
% with variation of perigee argument 
%
% Input: initial_orbit -> structure with initial orbital elements
%        final_orbit   -> structure with final orbital elements
%        f             -> low thrust acceleration
%        constants     ->

% Reference: Simplified Analysis of Low Thrust Orbital Manuevers, Pollard

% Marilena Di Carlo, 2015
% marilena.di-carlo@strath.ac.uk

% Assumption: continuous thrust all over the orbit (ref to Pollard)
alpha = pi/2;

% Variation of perigee argument
Delta_omega = final_orbit.omega - initial_orbit.omega;

% Equation 12 Pollard
DeltaV = 2 * alpha * Delta_omega / ...
    (sign(Delta_omega) * sqrt(initial_orbit.a / constants.mu_dim) * sqrt(1-initial_orbit.e^2) / initial_orbit.e * (3 * alpha - cos(alpha) * sin(alpha)) );

% Time of flight
ToF = DeltaV / f;