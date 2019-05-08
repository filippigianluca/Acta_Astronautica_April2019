function [Delta_omega, Delta_V] = omega0ToF_2_omegafDeltaV(initial_orbit, ToF, engine, constants)


% Input: initial_orbit
%        ToF
%        engine
%        constants
% Output: Delta_omega
%         Delta_V

% References: Petropolous, Ruggiero, Pollard Burt
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

%% Petropolous in plane

if initial_orbit.e ~= 0
    
    options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12);
    beta = 0;
    mode = 'optimum_omega';
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 ToF], ...
        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_omega_max_var_OOP = Kep_Elements(end,5) - initial_orbit.omega;
    Delta_V_max_var_OOP = Time(end) * engine.f;
    
    %% Ruggiero in plane out of plane
    
    options_ODE = odeset('AbsTol',1e-8,'RelTol',1e-8);
    beta = 0;
    mode = 'optimum_omega_Ruggiero';
    [Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 ToF], ...
        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
        options_ODE);
    
    Delta_omega_max_var_IP = Kep_Elements(end,5) - initial_orbit.omega;
    Delta_V_max_var_IP = Time(end) * engine.f;
    
    %% Burt aC
    
    
    Delta_omega_Burt_aC = 2/pi * sqrt(initial_orbit.a / constants.mu_dim) * ...
        (2 - initial_orbit.e^2) / initial_orbit.e * engine.f * ToF;
    
    
    Delta_V_Burt_aC = abs(Delta_omega_Burt_aC) / ...
        (2 / pi * sqrt(initial_orbit.a / constants.mu_dim) * (2 - initial_orbit.e^2) / initial_orbit.e);
    
    %% Burt aR
    
    Delta_omega_Burt_aR = sqrt(initial_orbit.a * (1 - initial_orbit.e^2) / constants.mu_dim) * ...
        engine.f * ToF;
    
    Delta_V_Burt_aR = abs(Delta_omega_Burt_aR) / ...
        sqrt(initial_orbit.a * (1 - initial_orbit.e^2) / constants.mu_dim);
    
    %% Pollard
    
    Delta_omega_Pollard = 3 * sqrt(initial_orbit.a * (1 - initial_orbit.e^2) / constants.mu_dim) * ...
        engine.f * ToF / (2 * initial_orbit.e);
    
    alpha = pi/2;
    Delta_V_Pollard = 2 * alpha * Delta_omega_Pollard / ...
        (sign(Delta_omega_Pollard) * sqrt(initial_orbit.a / constants.mu_dim) * sqrt(1-initial_orbit.e^2) / initial_orbit.e * (3 * alpha - cos(alpha) * sin(alpha)) );
    
    
    Delta_omega.max_var_IP  = Delta_omega_max_var_IP;
    Delta_omega.max_var_OOP = Delta_omega_max_var_OOP;
    Delta_omega.Burt_aC     = Delta_omega_Burt_aC;
    Delta_omega.Burt_aR     = Delta_omega_Burt_aR;
    Delta_omega.Pollard     = Delta_omega_Pollard;
    
    Delta_V.max_var_IP  = Delta_V_max_var_IP;
    Delta_V.max_var_OOP = Delta_V_max_var_OOP;
    Delta_V.Burt_aC     = Delta_V_Burt_aC;
    Delta_V.Burt_aR     = Delta_V_Burt_aR;
    Delta_V.Pollard     = Delta_V_Pollard;
    
else
    
    Delta_omega = [];
    Delta_V = [];
    
end
%%

% [max_omega_final, strategy] = max([max_omega_Burt_aC, max_omega_Burt_aR, max_omega_Pollard]);
% [min_omega_final, strategy] = min([min_omega_Burt_aC, min_omega_Burt_aR, min_omega_Pollard]);
%
% switch strategy
%     case 1
%         strategy_name = 'Burt aC';
%     case 2
%         strategy_name = 'Burt aR';
%     case 3
%         strategy_name = 'Pollard';
% end