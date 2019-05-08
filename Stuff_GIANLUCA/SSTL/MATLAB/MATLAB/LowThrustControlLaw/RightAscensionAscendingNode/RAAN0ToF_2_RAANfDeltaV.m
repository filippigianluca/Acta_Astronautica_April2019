function [Delta_RAAN, Delta_V] = RAAN0ToF_2_RAANfDeltaV(initial_orbit, ToF, engine, constants)


% Input: initial_orbit
%        ToF
%        engine
%        constants
% Output: Delta_RAAN
%         Delta_V

% References: Ruggiero, KEchichian
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

%% Ruggiero

if initial_orbit.e == 0

    Delta_RAAN_max_var =  2 * engine.f / pi * sqrt(initial_orbit.a / constants.mu_dim) * ...
        ToF / initial_orbit.i;
    
    Delta_V_max_var = pi/2 * sqrt(constants.mu_dim / initial_orbit.a) * sin(initial_orbit.i);


else 

    mode = 'Omega';
    element_target = 'Omega';
    options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.Omega));
    tic
    beta = 0;
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 ToF], ...
        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_RAAN_max_var = Kep_Elements(end,4) - Kep_Elements(1,4);
    Delta_V_max_var = Time(end) * engine.f;
end

%% Kechichian
if initial_orbit.e == 0

    Delta_RAAN_Kechichian =  2/(pi * sin(initial_orbit.i)) * ( atan( engine.f * ToF / v0) );

    Delta_V_Kechichian = sqrt(constants.mu_dim / initial_orbit.a) * ...
        sqrt(2 * (1 - cos(pi/2 * sin(initial_orbit.i) * abs(Delta_RAAN_Kechichian) )));

end

%% 

Delta_RAAN.max_var    = Delta_RAAN_max_var;
Delta_V.max_var       = Delta_V_max_var;

if initial_orbit.e == 0
    Delta_RAAN.Kechichian = Delta_RAAN_Kechichian;  
    Delta_V.Kechichian = Delta_V_Kechichian;
end

%%
% if initial_orbit.e 
%     
%     max_RAAN_final = max_RAAN_final_opt;
%     min_RAAN_final = min_RAAN_final_opt;
%     
%     strategy_name = 'Optimal';
%     
%     
% else
%     
%     [max_RAAN_final, strategy] = max([max_RAAN_final_opt, max_i_final_Kechichian]);
%     [min_RAAN_final, strategy] = min([min_RAAN_final_opt, min_i_final_Kechichian]);
%     
%     switch strategy
%         case 1
%             strategy_name = 'Optimal';
%         case 2
%             strategy_name = 'Kechichian';
%     end
%     
% end



