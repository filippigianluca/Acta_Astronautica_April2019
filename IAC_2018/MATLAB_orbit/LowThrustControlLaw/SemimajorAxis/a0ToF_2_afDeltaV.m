function [a_final_increase, a_final_decrease, DeltaV_increase, DeltaV_decrease] =...
    a0ToF_2_afDeltaV(initial_orbit, ...
    ToF, engine, constants)

% Input: initial_orbit
%        ToF
%        engine
%        constants
% Output: a_final_increase
%         a_final_decrease
%         DeltaV_increase
%         DeltaV_decrease

% References: Petropolous, Ruggiero, Edelbaum, Kechichian, Burt
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

%% Max instantaneous variation

% -------------------------------------------------------------------------
% Increase of semimajor axis
% -------------------------------------------------------------------------
mode = 'optimum_a';
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12);
beta = 0;
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 ToF], ...
                        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);

                    
DeltaV_increase.max_var = Time(end) * abs(engine.f);  
a_final_increase.max_var = Kep_Elements(end,1);


% -------------------------------------------------------------------------
% Decrease of semimajor axis
% -------------------------------------------------------------------------
mode = 'optimum_a';
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12);
engine.f = -engine.f;
engine.T = -engine.T;
beta = 0;
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 ToF], ...
                        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);

                    
DeltaV_decrease.max_var = Time(end) * abs(engine.f);  
a_final_decrease.max_var = Kep_Elements(end,1);



engine.f = -engine.f;
engine.T = -engine.T;
%%%%%%% APPROXIMATION:

[K,E] = ellipke(4 * initial_orbit.e / (1 + initial_orbit.e)^2);

fe = 2 / (1 - initial_orbit.e) * E + 2 / (1 + initial_orbit.e) * K;

% -------------------------------------------------------------------------
% Increase of semimajor axis
% -------------------------------------------------------------------------
% Final semimajor axis
a_final_increase.max_var_approx = ( 1 / initial_orbit.a + sign(engine.f) * ToF.^2 * engine.f^2 * ...
    (1 - initial_orbit.e^2)^2 * fe^2 / (4 * pi^2 * constants.mu_dim) + ...
    -  engine.f * (1 - initial_orbit.e^2) * fe * ToF / ...
    ( pi  * sqrt(constants.mu_dim * initial_orbit.a)) ).^(-1);

% Cost of the transfer
[K,E] = ellipke(4 * initial_orbit.e / (1 + initial_orbit.e)^2);
fe    = 2 / (1 - initial_orbit.e) * E + 2 / (1 + initial_orbit.e) * K;

DeltaV_increase.max_var_approx = 2 * pi * ...
    abs(sqrt(constants.mu_dim / initial_orbit.a) - sqrt(constants.mu_dim / a_final_increase.max_var_approx)) / ...
    (1 - initial_orbit.e^2) / fe;


% -------------------------------------------------------------------------
% Decrease of semimajor axis
% -------------------------------------------------------------------------
a_final_decrease.max_var_approx = ( 1 / initial_orbit.a + sign(engine.f) * ToF.^2 * engine.f^2 * ...
    (1 - initial_orbit.e^2)^2 * fe^2 / (4 * pi^2 * constants.mu_dim) + ...
    + engine.f * (1 - initial_orbit.e^2) * fe * ToF / ...
    ( pi  * sqrt(constants.mu_dim * initial_orbit.a)) ).^(-1);

% Cost of the transfer
[K,E] = ellipke(4 * initial_orbit.e / (1 + initial_orbit.e)^2);
fe    = 2 / (1 - initial_orbit.e) * E + 2 / (1 + initial_orbit.e) * K;

DeltaV_decrease.max_var_approx = 2 * pi * ...
    abs(sqrt(constants.mu_dim / initial_orbit.a) - sqrt(constants.mu_dim / a_final_decrease.max_var_approx)) / ...
    (1 - initial_orbit.e^2) / fe;


%% Edelbaum

if initial_orbit.e == 0

    v_initial = sqrt(constants.mu_dim / initial_orbit.a);

    % Increase of semimajor axis
    v_final_min = sqrt(v_initial^2 + engine.f^2 * ToF^2 - 2 * engine.f * ToF * v_initial);
    % Decrease of semimajor axis
    v_final_max = sqrt(v_initial^2 + engine.f^2 * ToF^2 + 2 * engine.f * ToF * v_initial);
    
    a_Edelbaum_increase = constants.mu_dim / v_final_min^2;
    a_Edelbaum_decrease = constants.mu_dim / v_final_max^2;
    
    DeltaV_Edelbaum_increase = abs(v_final_max - v_initial);
    DeltaV_Edelbaum_decrease = abs(v_final_min - v_initial);

end


%% Burt - check this equation

a_Burt_increase = initial_orbit.a * (1 - sqrt(initial_orbit.a / constants.mu_dim) * ...
    (2+initial_orbit.e^2 )* engine.f * ToF / (2 * sqrt(1-initial_orbit.e^2))).^(-2);

a_Burt_decrease = initial_orbit.a * (1 + sqrt(initial_orbit.a / constants.mu_dim) * ...
    (2+initial_orbit.e^2 )* engine.f * ToF / (2 * sqrt(1-initial_orbit.e^2))).^(-2);

DeltaV_Burt_increase = 2 * sqrt(constants.mu_dim) * sqrt(1-initial_orbit.e^2) * ...
    abs(1 / sqrt(initial_orbit.a) - 1 / sqrt(a_Burt_increase)) /...
           ( 2 + initial_orbit.e^2);
       
DeltaV_Burt_decrease = 2 * sqrt(constants.mu_dim) * sqrt(1-initial_orbit.e^2) * ...
    abs(1 / sqrt(initial_orbit.a) - 1 / sqrt(a_Burt_decrease)) /...
           ( 2 + initial_orbit.e^2);


%% Collect results in structure

% Increase of semimajor axis
% a_final_increase.max_var      = a_max_var_increase;
% DeltaV_increase.max_var = DeltaV_max_var_increase;

if initial_orbit.e == 0
    a_final_increase.Edelbaum = a_Edelbaum_increase;
    DeltaV_increase.Edelbaum = DeltaV_Edelbaum_increase;
end

a_final_increase.Burt     = a_Burt_increase;
DeltaV_increase.Burt      = DeltaV_Burt_increase;



% Decrease of semimajor axis
% a_final_decrease.max_var      = a_max_var_decrease;
% DeltaV_decrease.max_var = DeltaV_max_var_decrease;

if initial_orbit.e == 0
    a_final_decrease.Edelbaum = a_Edelbaum_decrease;
    DeltaV_decrease.Edelbaum = DeltaV_Edelbaum_decrease;
end

a_final_decrease.Burt     = a_Burt_decrease;
DeltaV_decrease.Burt      = DeltaV_Burt_decrease;





%% Compare resuts

% if initial_orbit.e == 0.1
%     [max_a_final, strategy] = max([a_optimal_max, a_Edelbaum_max, a_Burt_max]);
%     [min_a_final, strategy] = min([a_optimal_min, a_Edelbaum_min, a_Burt_min]);
%     
%     switch strategy
%         case 1
%             strategy_name = 'Petropolous';
%         case 2
%             strategy_name = 'Edelbaum';
%         case 3
%             strategy_name = 'Burt';
%     end
% else
%     [max_a_final, strategy] = max([a_optimal_max, a_Burt_max]);
%     [min_a_final, strategy] = min([a_optimal_min, a_Burt_min]);
%     
%     switch strategy
%         case 1
%             strategy_name = 'Petropolous';
%         case 2
%             strategy_name = 'Burt';
%     end
% end

