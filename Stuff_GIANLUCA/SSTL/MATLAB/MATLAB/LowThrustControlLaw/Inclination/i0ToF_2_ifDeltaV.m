function [Delta_i, Delta_V] = i0ToF_2_ifDeltaV(initial_orbit, ToF, engine, constants)

% Input: initial_orbit
%        ToF
%        engine
%        constants
% Output: Delta_i
%         Delta_V

% References: Petropolous, Edelbaum
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

%% Max variation

if initial_orbit.e == 0
    
    Delta_i.max_var =  2 * engine.f / pi * sqrt(initial_orbit.a / constants.mu_dim) * ToF;
    Delta_V.max_var = pi/2 * sqrt(constants.mu_dim / initial_orbit.a ) * Delta_i_max_var;
    
else
    
%     Delta_i.max_var_approx =  2 * engine.f / pi * sqrt(initial_orbit.a / constants.mu_dim) * ToF;
%     Delta_V.max_var_approx = pi/2 * sqrt(constants.mu_dim / initial_orbit.a ) * Delta_i_max_var;
%     
    
    mode = 'inclination';
    
    options_ODE = odeset('AbsTol',1e-10,'RelTol',1e-10);
    tic
    beta = 0;
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 ToF], ...
        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_V.max_var = Time(end) * abs(engine.f);
    Delta_i.max_var = Kep_Elements(end,3) - Kep_Elements(1,3);
    
end

%% Edelbaum

if initial_orbit.e == 0

    v0 = sqrt(constants.mu_dim / initial_orbit.a);
    Delta_i_Edelbaum =  2/pi * ( atan( engine.f * ToF / v0) );
    
    Delta_V_Edelbaum = v0 * sqrt(2 * (1 - cos(pi/2 * Delta_i_Edelbaum)));
end


%% 

% Delta_i.max_var  = Delta_i_max_var;
if initial_orbit.e == 0
    Delta_i.Edelbaum = Delta_i_Edelbaum;
end

% Delta_V.max_var  = Delta_V_max_var;
if initial_orbit.e == 0
    Delta_V.Edelbaum = Delta_V_Edelbaum;
end
%%
% if initial_orbit.e 
%     
%     max_i_final = max_i_final_opt;
%     min_i_final = min_i_final_opt;
%     
%     strategy_name = 'Optimal';
%     
%     
% else
%     
%     [max_i_final, strategy] = max([max_i_final_opt, max_i_final_Edelbaum]);
%     [min_i_final, strategy] = min([min_i_final_opt, min_i_final_Edelbaum]);
%     
%     switch strategy
%         case 1
%             strategy_name = 'Optimal';
%         case 2
%             strategy_name = 'Edelbaum';
%     end
%     
% end
