function [Delta_V, ToF] = RAAN0RAANf_2_DeltaVToF(initial_orbit, final_orbit, engine, constants)


% Input: initial_orbit
%        final_orbit
%        engine
%        constants
% Output: Delta_V
%         ToF

% References: Ruggiero, KEchichian
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

Delta_RAAN = abs(final_orbit.Omega - initial_orbit.Omega);


%% Ruggiero

if initial_orbit.e == 0 
    
    Delta_V.a0_e0_i0_max_var = pi/2 * sqrt(constants.mu_dim / initial_orbit.a) * sin(initial_orbit.i);
    ToF.a0_e0_i0_max_var     = Delta_V.a0_e0_i0_max_var;
    
    Delta_V.a0_e0_if_max_var = pi/2 * sqrt(constants.mu_dim / initial_orbit.a) * sin(final_orbit.i);
    ToF.a0_e0_if_max_var     = Delta_V.a0_e0_if_max_var;
    
    Delta_V.af_e0_i0_max_var = pi/2 * sqrt(constants.mu_dim / final_orbit.a) * sin(initial_orbit.i);
    ToF.af_e0_i0_max_var     = Delta_V.af_e0_i0_max_var;
    
    Delta_V.af_e0_if_max_var = pi/2 * sqrt(constants.mu_dim / final_orbit.a) * sin(final_orbit.i);
    ToF.af_e0_if_max_var     = Delta_V.af_e0_if_max_var;
    
    
else
    time_int = 365*86400;
    mode = 'Omega';
    element_target = 'Omega';
    options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.Omega));
    beta = 0;
    
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 time_int], ...
        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    Delta_V.a0_e0_i0_max_var = Time(end) * abs(engine.f);
    ToF.a0_e0_i0_max_var = Time(end);
    
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 time_int], ...
        [initial_orbit.a initial_orbit.e final_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    Delta_V.a0_e0_if_max_var = Time(end) * abs(engine.f);
    ToF.a0_e0_if_max_var = Time(end);
    
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 time_int], ...
        [final_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    Delta_V.af_e0_i0_max_var = Time(end) * abs(engine.f);
    ToF.af_e0_i0_max_var = Time(end);
    
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 time_int], ...
        [final_orbit.a initial_orbit.e final_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    Delta_V.af_e0_if_max_var = Time(end) * abs(engine.f);
    ToF.af_e0_if_max_var = Time(end);

end



if final_orbit.e == 0 
    
    Delta_V.a0_ef_i0_max_var = pi/2 * sqrt(constants.mu_dim / initial_orbit.a) * sin(initial_orbit.i);
    ToF.a0_ef_i0_max_var     = Delta_V.a0_ef_i0_max_var;
    
    Delta_V.a0_ef_if_max_var = pi/2 * sqrt(constants.mu_dim / initial_orbit.a) * sin(final_orbit.i);
    ToF.a0_ef_if_max_var     = Delta_V.a0_ef_if_max_var;
    
    Delta_V.af_ef_i0_max_var = pi/2 * sqrt(constants.mu_dim / final_orbit.a) * sin(initial_orbit.i);
    ToF.af_ef_i0_max_var     = Delta_V.af_ef_i0_max_var;
    
    Delta_V.af_ef_if_max_var = pi/2 * sqrt(constants.mu_dim / final_orbit.a) * sin(final_orbit.i);
    ToF.af_ef_if_max_var     = Delta_V.af_ef_if_max_var;
    
    
else
    time_int = 365*86400;
    mode = 'Omega';
    element_target = 'Omega';
    options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.Omega));
    beta = 0;
    
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 time_int], ...
        [initial_orbit.a final_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    Delta_V.a0_ef_i0_max_var = Time(end) * abs(engine.f);
    ToF.a0_ef_i0_max_var = Time(end);
    
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 time_int], ...
        [initial_orbit.a final_orbit.e final_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    Delta_V.a0_ef_if_max_var = Time(end) * abs(engine.f);
    ToF.a0_ef_if_max_var = Time(end);
    
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 time_int], ...
        [final_orbit.a final_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    Delta_V.af_ef_i0_max_var = Time(end) * abs(engine.f);
    ToF.af_ef_i0_max_var = Time(end);
    
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 time_int], ...
        [final_orbit.a final_orbit.e final_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    Delta_V.af_ef_if_max_var = Time(end) * abs(engine.f);
    ToF.af_ef_if_max_var = Time(end);
    
end


%% Kechichian
if initial_orbit.e == 0 


    Delta_V.a0_e0_i0_Kechichian = sqrt(constants.mu_dim / initial_orbit.a) * ...
        sqrt(2 * (1 - cos(pi/2 * sin(initial_orbit.i) * abs(Delta_RAAN) )));
    
    Delta_V.a0_e0_if_Kechichian = sqrt(constants.mu_dim / initial_orbit.a) * ...
        sqrt(2 * (1 - cos(pi/2 * sin(final_orbit.i) * abs(Delta_RAAN) )));
    
    Delta_V.af_e0_i0_Kechichian = sqrt(constants.mu_dim / final_orbit.a) * ...
        sqrt(2 * (1 - cos(pi/2 * sin(initial_orbit.i) * abs(Delta_RAAN) )));
    
    Delta_V.af_e0_if_Kechichian = sqrt(constants.mu_dim / final_orbit.a) * ...
        sqrt(2 * (1 - cos(pi/2 * sin(final_orbit.i) * abs(Delta_RAAN) )));
    

end


if final_orbit.e == 0 


    Delta_V.a0_ef_i0_Kechichian = sqrt(constants.mu_dim / initial_orbit.a) * ...
        sqrt(2 * (1 - cos(pi/2 * sin(initial_orbit.i) * abs(Delta_RAAN) )));
    
    Delta_V.af_e0_if_Kechichian = sqrt(constants.mu_dim / initial_orbit.a) * ...
        sqrt(2 * (1 - cos(pi/2 * sin(final_orbit.i) * abs(Delta_RAAN) )));
    
    Delta_V.af_ef_i0_Kechichian = sqrt(constants.mu_dim / final_orbit.a) * ...
        sqrt(2 * (1 - cos(pi/2 * sin(initial_orbit.i) * abs(Delta_RAAN) )));
    
    Delta_V.af_ef_if_Kechichian = sqrt(constants.mu_dim / final_orbit.a) * ...
        sqrt(2 * (1 - cos(pi/2 * sin(final_orbit.i) * abs(Delta_RAAN) )));
    

end





