function [Delta_e, Delta_V] = e0ToF_2_efDeltaV(initial_orbit, ToF, engine, constants)

% Input: initial_orbit
%        ToF
%        engine
%        constants
% Output: Delta_e
%         Delta_V

% References: Petropolous, Pollard, Burt
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

%% Petropolous

mode = 'optimum_eccentricity';
beta = 0;
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12);
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 ToF], ...
                        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);

e_max_var = Kep_Elements(end,2);

DeltaV_max_var = Time(end) * engine.f;

%% Pollard 

e_Pollard =  sin( asin(initial_orbit.e) + 1.5  * sqrt(initial_orbit.a / constants.mu_dim) * ...
    engine.f * ToF);

alpha = pi/2;
beta  = 0;
DeltaV_Pollard = sqrt(constants.mu_dim / initial_orbit.a) * ...
    (2 * alpha * abs(asin(initial_orbit.e) - asin(e_Pollard))) / ...
        (cos(abs(beta)) * (3 * alpha + cos(alpha) * sin(alpha)));


%% Burt

e_Burt = initial_orbit.e + sin(4/pi * sqrt(initial_orbit.a/constants.mu_dim) * engine.f * ToF );

DeltaV_Burt = pi/4 * sqrt(constants.mu_dim / initial_orbit.a) *...
    abs(asin(initial_orbit.e) - asin(e_Burt));    

%% 


Delta_e.max_var = e_max_var - initial_orbit.e;
Delta_e.Pollard = e_Pollard - initial_orbit.e;
Delta_e.Burt    = e_Burt    - initial_orbit.e;

Delta_V.max_var = DeltaV_max_var;
Delta_V.Pollard = DeltaV_Pollard;
Delta_V.Burt    = DeltaV_Burt;

%% 

% [max_e_final, strategy] = max([e_Pollard_max, e_Burt_max]);
% [min_e_final, strategy] = min([e_Pollard_min, e_Burt_min]);
% switch strategy
%     case 1
%         strategy_name = 'Pollard';
%     case 2
%         strategy_name = 'Burt';
% end


