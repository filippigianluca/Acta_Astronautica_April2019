function [Delta_V, ToF] = e0ef_2_DeltaVToF(initial_orbit, final_orbit, engine, constants)

% Input: initial_orbit
%        final_orbit
%        ToF
%        engine
%        constants
% Output: Delta_V
%         ToF

% References: Petropolous, Pollard, Burt
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

%% Petropolous



mode = 'optimum_eccentricity';
element_target = 'eccentricity';
beta = 0;
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.e));
f = sign(final_orbit.e - initial_orbit.e) * engine.f;
% -------------------------------------------------------------------------
% Variation of eccentricity with other orbital elements equal to their
% initial value
% -------------------------------------------------------------------------
if initial_orbit.e == 0
    e = 1e-10;
else
    e = initial_orbit.e;
end
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 365*86400], ...
                        [initial_orbit.a e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);

Delta_V.a0_max_var = Time(end) * engine.f;
ToF.a0_max_var = Time(end);


% -------------------------------------------------------------------------
% Variation of eccentricity with other orbital elements equal to their
% final value
% -------------------------------------------------------------------------
if initial_orbit.e == 0
    e = 1e-10;
else
    e = initial_orbit.e;
end
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x, f,beta,constants, mode),[0 365*86400], ...
                        [final_orbit.a e final_orbit.i final_orbit.Omega final_orbit.omega final_orbit.E], ...
                        options_ODE);

Delta_V.af_max_var = Time(end) * engine.f;
ToF.af_max_var = Time(end);

%% Pollard 

alpha = pi/2;
beta  = 0;

% -------------------------------------------------------------------------
% Variation of eccentricity with other orbital elements equal to their
% initial value
% -------------------------------------------------------------------------
Delta_V.a0_Pollard = sqrt(constants.mu_dim / initial_orbit.a) * ...
    (2 * alpha * abs(asin(initial_orbit.e) - asin(final_orbit.e))) / ...
        (cos(abs(beta)) * (3 * alpha + cos(alpha) * sin(alpha)));

ToF.a0_Pollard = Delta_V.a0_Pollard / engine.f;


% -------------------------------------------------------------------------
% Variation of eccentricity with other orbital elements equal to their
% final value
% -------------------------------------------------------------------------
Delta_V.af_Pollard = sqrt(constants.mu_dim / final_orbit.a) * ...
    (2 * alpha * abs(asin(initial_orbit.e) - asin(final_orbit.e))) / ...
        (cos(abs(beta)) * (3 * alpha + cos(alpha) * sin(alpha)));

ToF.af_Pollard = Delta_V.af_Pollard / engine.f;

%% Burt


% -------------------------------------------------------------------------
% Variation of eccentricity with other orbital elements equal to their
% initial value
% -------------------------------------------------------------------------
Delta_V.a0_Burt = pi/4 * sqrt(constants.mu_dim / initial_orbit.a) *...
    abs(asin(initial_orbit.e) - asin(final_orbit.e));    

ToF.a0_Burt = Delta_V.a0_Burt / engine.f;



% -------------------------------------------------------------------------
% Variation of eccentricity with other orbital elements equal to their
% final value
% -------------------------------------------------------------------------
Delta_V.af_Burt = pi/4 * sqrt(constants.mu_dim / final_orbit.a) *...
    abs(asin(initial_orbit.e) - asin(final_orbit.e));    

ToF.af_Burt = Delta_V.af_Burt / engine.f;





