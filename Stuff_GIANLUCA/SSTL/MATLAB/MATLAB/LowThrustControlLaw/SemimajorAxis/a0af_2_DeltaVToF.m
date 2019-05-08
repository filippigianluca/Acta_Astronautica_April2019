function [Delta_V ToF] = a0af_2_DeltaVToF(initial_orbit, ...
    final_orbit, engine, constants)

% Input: initial_orbit
%        final_orbit
%        engine
%        constants
% Output: Delta_V
%         ToF

% References: Petropolous, Ruggiero, Edelbaum, Kechichian, Burt
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

%% Max instantaneous variation

if final_orbit.a > initial_orbit.a
elseif final_orbit.a < initial_orbit.a
    engine.f = -engine.f;
    engine.T = -engine.T;
end

mode = 'optimum_a';

element_target = 'a';

time_int = 365*86400;
% -------------------------------------------------------------------------
% Variation of semimajor axis with other orbital elements equal to their
% initial value
% -------------------------------------------------------------------------
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.a));
beta = 0;
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 time_int], ...
                        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);

                    
Delta_V.e0_max_var = Time(end) * abs(engine.f);  
ToF.e0_max_var = Time(end);


% -------------------------------------------------------------------------
% Variation of semimajor axis with other orbital elements equal to their
% final value
% -------------------------------------------------------------------------
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.a));
beta = 0;
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 time_int], ...
                        [initial_orbit.a final_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);

                    
Delta_V.ef_max_var = Time(end) * abs(engine.f);  
ToF.ef_max_var = Time(end);


%%%%%%% APPROXIMATION:

% -------------------------------------------------------------------------
% Variation of semimajor axis with other orbital elements equal to their
% initial value
% -------------------------------------------------------------------------
[K,E] = ellipke(4 * initial_orbit.e / (1 + initial_orbit.e)^2);
fe    = 2 / (1 - initial_orbit.e) * E + 2 / (1 + initial_orbit.e) * K;

Delta_V.e0_max_var_approx = 2 * pi * ...
    abs(sqrt(constants.mu_dim / initial_orbit.a) - sqrt(constants.mu_dim / final_orbit.a)) / ...
    (1 - initial_orbit.e^2) / fe;

ToF.e0_max_var_approx = Delta_V.e0_max_var / engine.f;


% -------------------------------------------------------------------------
% Variation of semimajor axis with other orbital elements equal to their
% final value
% -------------------------------------------------------------------------
% Cost of the transfer
[K,E] = ellipke(4 * final_orbit.e / (1 + final_orbit.e)^2);
fe    = 2 / (1 - final_orbit.e) * E + 2 / (1 + final_orbit.e) * K;

Delta_V.ef_max_var_approx = 2 * pi * ...
    abs(sqrt(constants.mu_dim / initial_orbit.a) - sqrt(constants.mu_dim / final_orbit.a)) / ...
    (1 - final_orbit.e^2) / fe;

ToF.ef_max_var_approx = Delta_V.ef_max_var / engine.f;




%% Edelbaum

% -------------------------------------------------------------------------
% Variation of semimajor axis with other orbital elements equal to their
% initial value
% -------------------------------------------------------------------------
if initial_orbit.e == 0

    v_initial = sqrt(constants.mu_dim / initial_orbit.a);
    
    v_final = sqrt(constants.mu_dim / final_orbit.a);

    Delta_V.e0_Edelbaum = abs(v_final - v_initial);
    ToF.e0_Edelbaum = Delta_V.e0_Edelbaum / engine.f;

end


% -------------------------------------------------------------------------
% Variation of semimajor axis with other orbital elements equal to their
% final value
% -------------------------------------------------------------------------
if final_orbit.e == 0

    v_initial = sqrt(constants.mu_dim / initial_orbit.a);
    
    v_final = sqrt(constants.mu_dim / final_orbit.a);

    Delta_V.ef_Edelbaum = abs(v_final - v_initial);
    ToF.ef_Edelbaum = Delta_V.ef_Edelbaum / engine.f;

end


%% Burt - check this equation

% -------------------------------------------------------------------------
% Variation of semimajor axis with other orbital elements equal to their
% initial value
% -------------------------------------------------------------------------

Delta_V.e0_Burt = 2 * sqrt(constants.mu_dim) * sqrt(1-initial_orbit.e^2) * ...
    abs(1 / sqrt(initial_orbit.a) - 1 / sqrt(final_orbit.a)) /...
           ( 2 + initial_orbit.e^2);
       
ToF.e0_Burt = Delta_V.e0_Burt / engine.f;


% -------------------------------------------------------------------------
% Variation of semimajor axis with other orbital elements equal to their
% final value
% -------------------------------------------------------------------------

Delta_V.ef_Burt = 2 * sqrt(constants.mu_dim) * sqrt(1-final_orbit.e^2) * ...
    abs(1 / sqrt(initial_orbit.a) - 1 / sqrt(final_orbit.a)) /...
           ( 2 + final_orbit.e^2);
       
ToF.ef_Burt = Delta_V.ef_Burt / engine.f;
       






