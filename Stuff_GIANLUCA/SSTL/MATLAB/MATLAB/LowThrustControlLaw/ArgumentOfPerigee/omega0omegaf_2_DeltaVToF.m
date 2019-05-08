function [Delta_V, ToF ] = omega0omegaf_2_DeltaVToF(initial_orbit, final_orbit, engine, constants)


% Input: initial_orbit
%        final_orbit
%        engine
%        constants
% Output: Delta_V
%         ToF

% References: Petropolous, Ruggiero, Pollard Burt
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

Delta_omega = abs(final_orbit.omega - initial_orbit.omega);

if initial_orbit.e ~= 0
    
    %% Petropolous in plane
    
    mode = 'optimum_omega';
    element_target = 'omega';
    options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,...
        'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.omega));
    beta = 0;
    
    f = sign(final_orbit.omega - initial_orbit.omega) * engine.f;

    % -------------------------------------------------------------------------
    % Variation of omega with semimajor axis and eccentricity equal to their
    % initial values
    % -------------------------------------------------------------------------
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 365*86400], ...
        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_V.a0_e0_max_var_OOP = Time(end) * engine.f;
    ToF.a0_e0_max_var_OOP = Time(end);
    
    
    % -------------------------------------------------------------------------
    % Variation of omega with eccentricity equal to initial value and
    % semimajor axis equal to final value
    % initial values
    % -------------------------------------------------------------------------
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 365*86400], ...
        [final_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_V.af_e0_max_var_OOP = Time(end) * engine.f;
    ToF.af_e0_max_var_OOP = Time(end);
    
    %% Ruggiero in plane out of plane
    
    options_ODE = odeset('AbsTol',1e-8,'RelTol',1e-8, ...
        'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.omega));
    beta = 0;
    mode = 'optimum_omega_Ruggiero';
    
    % -------------------------------------------------------------------------
    % Variation of omega with semimajor axis and eccentricity equal to their
    % initial values
    % -------------------------------------------------------------------------
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 365*86400], ...
        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_V.a0_e0_max_var_IP = Time(end) * engine.f;
    
    
    % -------------------------------------------------------------------------
    % Variation of omega with eccentricity equal to initial value and
    % semimajor axis equal to final value
    % initial values
    % -------------------------------------------------------------------------
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 365*86400], ...
        [final_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_V.af_e0_max_var_IP = Time(end) * engine.f;
    
    %% Burt aC
 
    Delta_V.a0_e0_Burt_aC = Delta_omega / ...
        (2 / pi * sqrt(initial_orbit.a / constants.mu_dim) * (2 - initial_orbit.e^2) / initial_orbit.e);
    
    Delta_V.af_e0_Burt_aC = Delta_omega / ...
        (2 / pi * sqrt(final_orbit.a / constants.mu_dim) * (2 - initial_orbit.e^2) / initial_orbit.e);

    
    %% Burt aR
    
    Delta_V.a0_e0_Burt_aR = Delta_omega / ...
        sqrt(initial_orbit.a * (1 - initial_orbit.e^2) / constants.mu_dim);
    
    Delta_V.af_e0_Burt_aR = Delta_omega / ...
        sqrt(final_orbit.a * (1 - initial_orbit.e^2) / constants.mu_dim);
    
    %% Pollard
    
 
    alpha = pi/2;
    Delta_V.a0_e0_Pollard = 2 * alpha * Delta_omega / ...
        (sign(Delta_omega) * sqrt(initial_orbit.a / constants.mu_dim) * sqrt(1-initial_orbit.e^2) / initial_orbit.e * (3 * alpha - cos(alpha) * sin(alpha)) );
    
    alpha = pi/2;
    Delta_V.af_e0_Pollard = 2 * alpha * Delta_omega / ...
        (sign(Delta_omega) * sqrt(final_orbit.a / constants.mu_dim) * sqrt(1-initial_orbit.e^2) / initial_orbit.e * (3 * alpha - cos(alpha) * sin(alpha)) );
    
    
        
end






if final_orbit.e ~= 0
    
    %% Petropolous in plane
    
    mode = 'optimum_omega';
    element_target = 'omega';
    options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,...
        'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.omega));
    beta = 0;
    
    f = sign(final_orbit.omega - initial_orbit.omega) * engine.f;

    % -------------------------------------------------------------------------
    % Variation of omega with eccentricity equal to final value and
    % semimajor axis equal to initial value
    % -------------------------------------------------------------------------
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 365*86400], ...
        [initial_orbit.a final_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_V.a0_ef_max_var_OOP = Time(end) * engine.f;
    ToF.a0_ef_max_var_OOP = Time(end);
    
    
    % -------------------------------------------------------------------------
    % Variation of omega with eccentricity equal to final values
    % -------------------------------------------------------------------------
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 365*86400], ...
        [final_orbit.a final_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_V.af_ef_max_var_OOP = Time(end) * engine.f;
    ToF.af_ef_max_var_OOP = Time(end);
    
    %% Ruggiero in plane out of plane
    
    options_ODE = odeset('AbsTol',1e-8,'RelTol',1e-8, ...
        'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.omega));
    beta = 0;
    mode = 'optimum_omega_Ruggiero';
    
    % -------------------------------------------------------------------------
    % Variation of omega with semimajor axis equal to initial value and
    % eccentricity equal to final value
    % -------------------------------------------------------------------------
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 365*86400], ...
        [initial_orbit.a final_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_V.a0_ef_max_var_IP = Time(end) * engine.f;
    
    
    % -------------------------------------------------------------------------
    % Variation of omega with eccentricity and semimajor axis equal to
    % final value
    % -------------------------------------------------------------------------
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 365*86400], ...
        [final_orbit.a final_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_V.af_ef_max_var_IP = Time(end) * engine.f;
    
    %% Burt aC
 
    Delta_V.a0_ef_Burt_aC = Delta_omega / ...
        (2 / pi * sqrt(initial_orbit.a / constants.mu_dim) * (2 - final_orbit.e^2) / final_orbit.e);
    
    Delta_V.af_ef_Burt_aC = Delta_omega / ...
        (2 / pi * sqrt(final_orbit.a / constants.mu_dim) * (2 - final_orbit.e^2) / final_orbit.e);

    
    %% Burt aR
    
    Delta_V.a0_ef_Burt_aR = Delta_omega / ...
        sqrt(initial_orbit.a * (1 - final_orbit.e^2) / constants.mu_dim);
    
    Delta_V.af_ef_Burt_aR = Delta_omega / ...
        sqrt(final_orbit.a * (1 - final_orbit.e^2) / constants.mu_dim);
    
    %% Pollard
    
 
    alpha = pi/2;
    Delta_V.a0_ef_Pollard = 2 * alpha * Delta_omega / ...
        (sign(Delta_omega) * sqrt(initial_orbit.a / constants.mu_dim) * sqrt(1-final_orbit.e^2) / final_orbit.e * (3 * alpha - cos(alpha) * sin(alpha)) );
    
    alpha = pi/2;
    Delta_V.af_ef_Pollard = 2 * alpha * Delta_omega / ...
        (sign(Delta_omega) * sqrt(final_orbit.a / constants.mu_dim) * sqrt(1-final_orbit.e^2) / final_orbit.e * (3 * alpha - cos(alpha) * sin(alpha)) );
    
    
        
end