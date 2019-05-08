function [Delta_V, ToF] = i0if_2_DeltaVToF(initial_orbit, final_orbit, engine, constants)

% Input: initial_orbit
%        final_orbit
%        engine
%        constants
% Output: DeltaV
%         ToF

% References: Petropolous, Edelbaum
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk


Delta_i =  abs(final_orbit.i - initial_orbit.i);

%% Max variation

% -------------------------------------------------------------------------
% Variation of inclination with other orbital elements equal to their
% initial value
% -------------------------------------------------------------------------
if initial_orbit.e == 0
    Delta_V.a0_e0_max_var = pi/2 * sqrt(constants.mu_dim / initial_orbit.a ) * Delta_i;
    ToF.a0_e0_max_var = Delta_V.a0_e0_max_var / engine.f;
else
    mode = 'inclination';
    element_target = 'inclination';
    options_ODE = odeset('AbsTol',1e-10,'RelTol',1e-10,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.i));
    tic
    beta = 0;
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 365*86400], ...
        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_V.a0_e0_max_var = Time(end) * abs(engine.f);
    ToF.a0_e0_max_var     = Time(end);
end


% -------------------------------------------------------------------------
% Variation of inclination with semimajor axis equal to final value and
% eccentricity equal to initial value
% -------------------------------------------------------------------------
if initial_orbit.e == 0
    Delta_V.af_e0_max_var = pi/2 * sqrt(constants.mu_dim / final_orbit.a ) * Delta_i;
    ToF.af_e0_max_var = Delta_V.af_e0_max_var / engine.f;
else
    mode = 'inclination';
    element_target = 'inclination';
    options_ODE = odeset('AbsTol',1e-10,'RelTol',1e-10,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.i));
    tic
    beta = 0;
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 365*86400], ...
        [final_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_V.af_e0_max_var = Time(end) * abs(engine.f);
    ToF.af_e0_max_var     = Time(end);
end


% -------------------------------------------------------------------------
% Variation of inclination with other orbital elements equal to their
% final value
% -------------------------------------------------------------------------
if final_orbit.e == 0
    Delta_V.af_ef_max_var = pi/2 * sqrt(constants.mu_dim / final_orbit.a ) * Delta_i;
    ToF.af_ef_max_var = Delta_V.af_ef_max_var / engine.f;

else

    mode = 'inclination';
    element_target = 'inclination';
    options_ODE = odeset('AbsTol',1e-10,'RelTol',1e-10,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.i));
    tic
    beta = 0;
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 365*86400], ...
        [final_orbit.a final_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_V.af_ef_max_var = Time(end) * abs(engine.f);
    ToF.af_ef_max_var     = Time(end);
    
end



% -------------------------------------------------------------------------
% Variation of inclination with eccentricity equal to final value and
% semimajor axis equal to initial value
% -------------------------------------------------------------------------
if final_orbit.e == 0
    Delta_V.a0_ef_max_var = pi/2 * sqrt(constants.mu_dim / initial_orbit.a ) * Delta_i;
    ToF.a0_ef_max_var = Delta_V.a0_ef_max_var / engine.f;

else

    mode = 'inclination';
    element_target = 'inclination';
    options_ODE = odeset('AbsTol',1e-10,'RelTol',1e-10,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.i));
    tic
    beta = 0;
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 365*86400], ...
        [initial_orbit.a final_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_V.a0_ef_max_var = Time(end) * abs(engine.f);
    ToF.a0_ef_max_var     = Time(end);
    
end

%% Edelbaum

if initial_orbit.e == 0
    
    % -------------------------------------------------------------------------
    % Variation of inclination with other orbital elements equal to their
    % initial value
    % -------------------------------------------------------------------------
    
    v0 = sqrt(constants.mu_dim / initial_orbit.a);
    
    Delta_V.a0_Edelbaum = v0 * sqrt(2 * (1 - cos(pi/2 * Delta_i)));
    ToF.a0_Edelbaum     = Delta_V.a0_Edelbaum / engine.f;
    
    
    % -------------------------------------------------------------------------
    % Variation of inclination with other orbital elements equal to their
    % final value
    % -------------------------------------------------------------------------
    
    v0 = sqrt(constants.mu_dim / final_orbit.a);
    
    Delta_V.af_Edelbaum = v0 * sqrt(2 * (1 - cos(pi/2 * Delta_i)));
    ToF.af_Edelbaum     = Delta_V.af_Edelbaum / engine.f;
end



