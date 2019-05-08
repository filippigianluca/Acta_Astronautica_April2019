% =========================================================================
% Function for the computation of the objective function to be minimized
% and of the non-linear equality constraints
% =========================================================================
% Input: x          -> vector to optimise
%        parameters ->
%        options    ->
%        constants  ->

% Output: J         -> function to minimise (time spent with engine on)
%         C         -> non-linear constraints
%         Ceq       -> linear constraints

% The vector x is 1x14*arcs+1. 
% For each arcs 9 variable are defined: 
% - alpha (azimuth)
% - beta(elevation)
% - 6 equinoctial variables for the initial point of the propulsed arc 
% - 6 equinoctial variabes for the final point of the propulsed arc.

% They are organized in the x vector as follows:
% [alpha, beta, [a, P1, P2, Q1, Q2, L], [a, P1, P2, Q1, Q2, Delta L]]
% for each propulsed arc. 

% The last element in the vector x is the time of flight.

% Author: Marilena Di Carlo
% email: marilena.di-carlo@strath.ac.uk

function [J_out,C,Ceq,GCeq,ToF,Equin_sorted,Equin_sorted_ct,J] = ...
    cost_constraint_function_LT_2(x, parameters, options, constants)

% Parameters and variable for the transfer
arcs         = parameters.arcs;
eps_max      = parameters.eps_max;
departure_eq   = parameters.departure_eq;
% curr_tof     = parameters.curr_tof;
n            = parameters.n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Keplerian elements of arrival body at time 0
arrival_kep_t0 = parameters.arrival_kep_t0;

% Initial true anomaly
theta0 = arrival_kep_t0(6);

% Initial eccentric anomaly
cos_E0 = (cos(theta0) + arrival_kep_t0(2)) / ...
    (1 + arrival_kep_t0(2) * cos(theta0));

sin_E0 = (sin(theta0) * sqrt(1 - arrival_kep_t0(2)^2)) / ...
    (1 + arrival_kep_t0(2) * cos(theta0));

E0 = atan2(sin_E0, cos_E0);

% Initial mean anomaly
M0 = E0 - arrival_kep_t0(2) * sin(E0);

% Arrival time
ToF_opt = x(end);
ToF_opt * constants.TU;

% Solve Kepler problem
n_mean_mot = sqrt(constants.mu / arrival_kep_t0(1)^3);

% Orbital period arrival orbit
T = 2 * pi / n_mean_mot;

k = fix(ToF_opt/T);

% Final mean anomaly
M = M0 + n_mean_mot * ToF_opt - 2 * k * pi;

% Solve kepler problem to find final eccentric anomaly
E = kepler_E(arrival_kep_t0(2), M);

% Final true anomaly
cos_theta = (cos(E) - arrival_kep_t0(2)) /...
    (1 - arrival_kep_t0(2) * cos(E));

sin_theta = sin(E) * sqrt(1 - arrival_kep_t0(2)^2) / ...
        (1 - arrival_kep_t0(2) * cos(E));
    
thetaf = atan2(sin_theta, cos_theta);

% Final keplerian elements arrival body
arrival_kep_tf = arrival_kep_t0;
arrival_kep_tf(6) = thetaf;
thetaf(end);
% Equinoctial elements arrival body
arrival_eq = kep2eq(arrival_kep_tf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END NEW SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cost function
J = 0;

% Equality constraints
Ceq = [];



% Time of flight
ToF = 0;

% Initialize varaibles
Equin_coast = [];
Equin_thrust = [];




%% FIRST PROPULSED ARC

% Initial equinoctial element for forward propagation on the first low
% thrust arc
Equin_initial_prop_forw  = x(3:8);

% Final equinoctial elements of first forward propagation on the first low
% thrust arc
Equin_final_prop_forw = [x(9:13) x(8)+x(14)];

% True longitude variation over first low thrust arc
L_prop_forw  = linspace(Equin_initial_prop_forw(6),  Equin_final_prop_forw(6), n);


%% LAST PROPULSED ARC

% Initial equinoctial element for backward propagation on the first LT
% backward propagated arc (last LT arcs)
Equin_initial_prop_back  = [x(end-6:end-2) x(end-7)+x(end-1)];

% Final equinoctial elements for backward propagation on the first LT
% backward propagated arc (last LT arcs)
Equin_final_prop_back = x(end-12:end-7);

% Longitude values for propagation
L_prop_back  = linspace(Equin_initial_prop_back(6),  Equin_final_prop_back(6), n);

%% Initialise GCeq

% if options.gradient
    if (departure_eq(6) == Equin_initial_prop_forw(6)) && (arrival_eq(6) == Equin_initial_prop_back(6))
        rows_GCeq = 5 * (2 * arcs - 1);
    elseif ((departure_eq(6) == Equin_initial_prop_forw(6)) && (arrival_eq(6) ~= Equin_initial_prop_back(6))) || ...
            ((departure_eq(6) ~= Equin_initial_prop_forw(6)) && (arrival_eq(6) == Equin_initial_prop_back(6)))
        rows_GCeq = 5 * 2 * arcs;
    elseif ((departure_eq(6) ~= Equin_initial_prop_forw(6)) && (arrival_eq(6) ~= Equin_initial_prop_back(6)))
        rows_GCeq = 5 * (2*arcs + 1);
    end
    
    col_GCeq = length(x);
    
    GCeq = zeros(rows_GCeq, col_GCeq);
    
    delta = 1e-8;
    delta_a = 1e-8;
    adim_sma = 1;
% end
%% FIRST COAST ARC?

% Define an initial coast arc only if the first point of the vector x
% does not coincide with the departure point

if departure_eq(6) == Equin_initial_prop_forw(6)
    
    % If the longitude of the initial point of the first LT arc coincides with the
    % departure longitude, there is no need of using a coast arc and the
    % first equality contraint (coincidence of the two points) has to be
    % defined

    Ceq = [Ceq; departure_eq(1:5) - Equin_initial_prop_forw(1:5)'];
   
    GCeq(1,3) = -1;
    GCeq(2,4) = -1;
    GCeq(3,5) = -1;
    GCeq(4,6) = -1;
    GCeq(5,7) = -1;
    
    
    
else
    % If the initial point of the first propulsed arc does not coincide
    % with the departure point, use a coast arc to reach the first point
    
    
    % True longitude variation
    L_coast_forw = linspace(departure_eq(6), Equin_initial_prop_forw(6), n);
    
%     [Equin_coast_forw, t] = AnEquin_all_forward_tang_m(L_coast_forw, departure_eq, ...
%                                                        0, 0, ...
%                                                        0, 0, 0, ...
%                                                        0, 0, 0, ...
%                                                        constants.mu, options.geopotential, options.third_body,...
%                                                        constants.R_Earth, options.drag, 0);
    
    if strcmp(options.fun, 'AnEquin_forward_m')
        [Equin_coast_forw, t] = AnEquin_forward_m(L_coast_forw, departure_eq, 0, 0, 0, constants.mu);
    elseif strcmp(options.fun, 'AnEquin_all_forward_tang_m')
        [Equin_coast_forw, t] = AnEquin_all_forward_tang_m(L_coast_forw, departure_eq, 0, 0, ...
            0, 0, 0, 0, 0, 0, constants.mu, options.geopotential, options.third_body, constants.R_Earth, options.drag, 0);
    end

    Equin_coast = [Equin_coast, Equin_coast_forw];
    
    Equin_sorted{1} = Equin_coast_forw;
    Equin_sorted_ct(1) = 0;
    
    % Time of flight update
    ToF = ToF + t(end);
    
    % The final point of the initial coast arc has to coincide with the initial point
    % of the first low thrust arc

    Ceq = [Ceq; Equin_coast_forw(1:5,end) - Equin_initial_prop_forw(1:5)'];
    
    if options.gradient == 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Jacobian constraints
        x_var = x;
        x_var(8) = x_var(8) + delta;
        
        eps_GCeq = 0;
        alpha_GCeq = 0;
        beta_GCeq = 0;
        L_prop_GCeq = linspace(departure_eq(6), departure_eq(6)+x_var(8), n);
        start_eq_GCeq = departure_eq;
        check_eq_GCeq = Equin_initial_prop_forw(1:5)';
        check_eq_GCeq(1) = check_eq_GCeq(1);
        Ceq_var = jacobian_const_LT_OOP(eps_GCeq, alpha_GCeq, beta_GCeq, L_prop_GCeq, start_eq_GCeq, check_eq_GCeq, options, constants);
        Ceq_var(1) = Ceq_var(1);
        GCeq(1:5,8) = (Ceq_var - Ceq(1:5)) / delta;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        GCeq(1,3) = -1;
        GCeq(2,4) = -1;
        GCeq(3,5) = -1;
        GCeq(4,6) = -1;
        GCeq(5,7) = -1;
    end
    
    
end





%% LAST COAST ARC?


if arrival_eq(6) == Equin_initial_prop_back(6)
%        Ceq = [Ceq; (arrival_eq(1) - Equin_initial_prop_back(1)' ) / adim_sma;arrival_eq(2:5) - Equin_initial_prop_back(2:5)'];
       Ceq = [Ceq; arrival_eq(1:5) - Equin_initial_prop_back(1:5)'];
       
       GCeq(6,end-5) = -1;
       GCeq(7,end-4) = -1;
       GCeq(8,end-3) = -1;
       GCeq(9,end-2) = -1;
       GCeq(10,end-1) = -1;
else
    
    L_coast_back = linspace(arrival_eq(6), Equin_initial_prop_back(6), n);
    
%     [Equin_coast_back, t] = AnEquin_all_forward_tang_m(L_coast_back, arrival_eq, ...
%                                                         0, 0, ...
%                                                         0, 0, 0, ...
%                                                         0, 0, 0, ...
%                                                         constants.mu, options.geopotential, options.third_body, ...
%                                                         constants.R_Earth, options.drag, 0);


    if strcmp(options.fun, 'AnEquin_forward_m')
        [Equin_coast_back, t] = AnEquin_forward_m(L_coast_back, arrival_eq, 0, 0, 0, constants.mu);
    elseif strcmp(options.fun, 'AnEquin_all_forward_tang_m')
        [Equin_coast_back, t] = AnEquin_all_forward_tang_m(L_coast_back, arrival_eq, 0, 0, ...
            0, 0, 0, 0, 0, 0, constants.mu, options.geopotential, options.third_body, constants.R_Earth, options.drag, 0);
    end
    
    Equin_coast = [Equin_coast, Equin_coast_back];
    
    Equin_sorted{2*arcs+1} = Equin_coast_back(:,end:-1:1);
    Equin_sorted_ct(2*arcs+1) = 0;
    
    % Time of flight update (negative time for backward propagation)
    ToF = ToF - t(end);
    
    % Final point of the initial coast arc to coincide with the initial point
    % of the first backward propagated LT arc
    Ceq = [Ceq; Equin_coast_back(1:5,end) - Equin_initial_prop_back(1:5)'];
    
    if options.gradient == 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gradient constraints
        x_var = x;
        x_var(end) = x_var(end) + delta;
        
        %     [~,~,Ceq_var] = cost_constraint_function_OutOfPlane_2nd(x_var, parameters, options, constants);
        
        eps_GCeq = 0;
        alpha_GCeq = 0;
        beta_GCeq = 0;
        L_prop_GCeq = linspace(arrival_eq(6), x_var(end-6)+x_var(end), n);
        start_eq_GCeq = arrival_eq;
        check_eq_GCeq = Equin_initial_prop_back(1:5)';
        check_eq_GCeq(1) = check_eq_GCeq(1) ;
        Ceq_var = jacobian_const_LT_OOP(eps_GCeq,  alpha_GCeq, beta_GCeq, L_prop_GCeq, start_eq_GCeq, check_eq_GCeq, options, constants);
        Ceq_var(1) = Ceq_var(1);
        GCeq(6:10,end) = (Ceq_var - Ceq(6:10)) / delta;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        GCeq(6,length(x)-5) = -1;
        GCeq(7,length(x)-4) = -1;
        GCeq(8,length(x)-3) = -1;
        GCeq(9,length(x)-2) = -1;
        GCeq(10,length(x)-1) = -1;
    end
end

temp = 11;

%%

for i = 1 : arcs/2
%     if i == arcs/2
%         keyboard
%     end
    
    
    % ---------------------------------------------------------------------
    % Forward propagated low-thrust arc
    % ---------------------------------------------------------------------
    
    % Forward propagation from initial to final point of LT arc using
    % LT propulsion
%     [EquinLT_forward, t] = AnEquin_all_forward_tang_m(L_prop_forw, Equin_initial_prop_forw', ...
%                                                       0, 0, ...
%                                                       eps_max, x(1+(i-1)*9), 0, ...
%                                                       0, 0, 0, ...
%                                                       constants.mu, options.geopotential, options.third_body,...
%                                                       constants.R_Earth, options.drag, 0);

%     [EquinLT_forward, t] = AnEquin_forward_m(L_prop_forw, Equin_initial_prop_forw', eps_max, x(1+(i-1)*14), x(2+(i-1)*14), constants.mu);
    

    if strcmp(options.fun, 'AnEquin_forward_m')
        [EquinLT_forward, t] = AnEquin_forward_m(L_prop_forw, Equin_initial_prop_forw', eps_max, x(1+(i-1)*14), x(2+(i-1)*14), constants.mu);
    elseif strcmp(options.fun, 'AnEquin_all_forward_tang_m')
        [EquinLT_forward, t] = AnEquin_all_forward_tang_m(L_prop_forw, Equin_initial_prop_forw', 0, 0, ...
                                eps_max, x(1+(i-1)*14), x(2+(i-1)*14), 0, 0, 0, constants.mu, ...
                                options.geopotential, options.third_body, constants.R_Earth, options.drag, 0);
    end
    

    Equin_thrust = [Equin_thrust, EquinLT_forward];
    
    Equin_sorted{2+(i-1)*2} = EquinLT_forward;
    Equin_sorted_ct(2+(i-1)*2) = 1;
    
    % Cost function: time spent using propulsion
    J = J + t(end);
%     if t(end)<0
%         keyboard
%     end
%         
    
    % Time of flight update
    ToF = ToF + t(end);
    
    
    % Matching condition between final point of forward propagated propulsed
    % arc and required final point of forward propagated propulsed arc
    Ceq = [Ceq; EquinLT_forward(1:5,end) - Equin_final_prop_forw(1:5)'];
    
    if options.gradient == 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gradient constraints
    index_GCeq = [1 2 3 4 5 6 7 8 14] + (i-1)*14;
    eps_GCeq = eps_max;
    for i_GCeq = 1 : length(index_GCeq)
        x_var = x;
        if i_GCeq == 3 +  (i-1)*14
            delta_temp = delta_a;
        else
            delta_temp = delta;
        end
        x_var(index_GCeq(i_GCeq)) = x_var(index_GCeq(i_GCeq)) + delta_temp;
        
        %         [~,~,Ceq_var] = cost_constraint_function_OutOfPlane_2nd(x_var, parameters, options, constants);
        
        alpha_GCeq = x_var(1+(i-1)*14);
        beta_GCeq  = x_var(2+(i-1)*14);
        %         L_prop_GCeq = linspace(Equin_initial_prop_forw(6),  Equin_final_prop_forw(6), n);
        start_eq_GCeq = x_var(3+14*(i-1): 8+14*(i-1));
        start_eq_GCeq(1) = start_eq_GCeq(1);
        
        check_eq_GCeq = [x_var(9+14*(i-1): 13+14*(i-1))'; x_var(8+14*(i-1))+x_var(14+14*(i-1))];
        check_eq_GCeq(1) = check_eq_GCeq(1);
        
        L_prop_GCeq = linspace(start_eq_GCeq(6), check_eq_GCeq(6), n);
        
        Ceq_var = jacobian_const_LT_OOP(eps_GCeq, alpha_GCeq, beta_GCeq, L_prop_GCeq, start_eq_GCeq, check_eq_GCeq(1:5), options, constants);
        Ceq_var(1) = Ceq_var(1);
        
        GCeq(temp,index_GCeq(i_GCeq)) = (Ceq_var(1) - Ceq(temp)) / delta_a;
        GCeq(temp+1:temp+4,index_GCeq(i_GCeq)) = (Ceq_var(2:end) - Ceq(temp+1:temp+4)) / delta;
        
    end
    GCeq(temp,   9 + (i-1)*14) = -1;
    GCeq(temp+1, 10 + (i-1)*14) = -1;
    GCeq(temp+2, 11 + (i-1)*14) = -1;
    GCeq(temp+3, 12 + (i-1)*14) = -1;
    GCeq(temp+4, 13 + (i-1)*14) = -1;
    
    temp = temp + 5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
        
    % ---------------------------------------------------------------------
    % Forward coast arc from final point of low thurst arc to initial point
    % of next propulsed arc
    % ---------------------------------------------------------------------
    
    % Forward coast arc
    L_coast_forw = linspace(Equin_final_prop_forw(6), x(22+14*(i-1)), n);
    Equin_initial_coast_forw = Equin_final_prop_forw';
    
%     [Equin_coast_forward, t] = AnEquin_all_forward_tang_m(L_coast_forw, Equin_initial_coast_forw, ...
%                                                           0, 0, ...
%                                                           0, 0, 0, ...
%                                                           0, 0, 0, ...
%                                                           constants.mu, options.geopotential, options.third_body,...
%                                                           constants.R_Earth, options.drag, 0);

%     [Equin_coast_forward, t] = AnEquin_forward_m(L_coast_forw, Equin_initial_coast_forw, 0, 0, 0, constants.mu);

    if strcmp(options.fun, 'AnEquin_forward_m')
        [Equin_coast_forward, t] = AnEquin_forward_m(L_coast_forw, Equin_initial_coast_forw, 0, 0, 0, constants.mu);
    elseif strcmp(options.fun, 'AnEquin_all_forward_tang_m')
        [Equin_coast_forward, t] = AnEquin_all_forward_tang_m(L_coast_forw, Equin_initial_coast_forw, 0, 0, ...
                                0, 0, 0, 0, 0, 0, constants.mu, ...
                                options.geopotential, options.third_body, constants.R_Earth, options.drag, 0);
    end
    
    Equin_coast = [Equin_coast, Equin_coast_forward];
    
    Equin_sorted{3+(i-1)*2} = Equin_coast_forward;
    Equin_sorted_ct(3+(i-1)*2) = 0;
    
    % Time of flight update
    ToF = ToF + t(end);
    
    % Matching condition between final point of forward coast arc and
    % initial point of next propulsed arc

    Ceq = [Ceq; Equin_coast_forward(1:5,end) - x(17+14*(i-1):21+14*(i-1))'];
    
    if options.gradient == 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gradient constraints
    index_GCeq =  [9 10 11 12 13 14 22] + (i-1)*14;
    eps_GCeq = 0;
    for i_GCeq = 1 : length(index_GCeq)
        x_var = x;
        
        if i_GCeq == 9 +  (i-1)*14
            delta_temp = delta_a;
        else
            delta_temp = delta;
        end
        
        x_var(index_GCeq(i_GCeq)) = x_var(index_GCeq(i_GCeq)) + delta_temp;
        
%         [~,~,Ceq_var] = cost_constraint_function_OutOfPlane_2nd(x_var, parameters, options, constants);

        alpha_GCeq = 0;
        beta_GCeq  = 0;
        start_eq_GCeq = [x_var(9+14*(i-1): 13+14*(i-1)) x_var(8+14*(i-1))+x_var(14+14*(i-1))];
        start_eq_GCeq(1) = start_eq_GCeq(1) ;
        check_eq_GCeq = x_var(17 + 14 * (i-1) : 22+14*(i-1))';
        check_eq_GCeq(1) = check_eq_GCeq(1);
        
        L_prop_GCeq = linspace(start_eq_GCeq(6), check_eq_GCeq(6), n);

        Ceq_var = jacobian_const_LT_OOP(eps_GCeq, alpha_GCeq, beta_GCeq, L_prop_GCeq, start_eq_GCeq, check_eq_GCeq(1:5), options, constants);
        Ceq_var(1) = Ceq_var(1);
        GCeq(temp,index_GCeq(i_GCeq)) = (Ceq_var(1) - Ceq(temp)) / delta_a;
        GCeq(temp+1:temp+4,index_GCeq(i_GCeq)) = (Ceq_var(2:end) - Ceq(temp+1:temp+4)) / delta;
        
    end
    
    GCeq(temp,   17 + (i-1)*14) = -1;
    GCeq(temp+1, 18 + (i-1)*14) = -1;
    GCeq(temp+2, 19 + (i-1)*14) = -1;
    GCeq(temp+3, 20 + (i-1)*14) = -1;
    GCeq(temp+4, 21 + (i-1)*14) = -1;
    
    temp = temp + 5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    % ---------------------------------------------------------------------
    % Backward propagated low-thrust arc
    % ---------------------------------------------------------------------
    % Backward propagation from final to initial point of LT arc
%     [EquinLT_backward, t] = AnEquin_all_forward_tang_m(L_prop_back, Equin_initial_prop_back', ...
%                                                        0, 0, ...
%                                                        eps_max, x(end-8-9*(i-1)), 0, ...
%                                                        0, 0, 0, ...
%                                                        constants.mu, options.geopotential, options.third_body,...
%                                                        constants.R_Earth, options.drag, 0);

%     [EquinLT_backward, t] = AnEquin_forward_m(L_prop_back, Equin_initial_prop_back', eps_max, x(end-13-14*(i-1)), x(end-12-14*(i-1)), constants.mu);

    if strcmp(options.fun, 'AnEquin_forward_m')
        [EquinLT_backward, t] = AnEquin_forward_m(L_prop_back, Equin_initial_prop_back', eps_max, x(end-14-14*(i-1)), x(end-13-14*(i-1)), constants.mu);
    elseif strcmp(options.fun, 'AnEquin_all_forward_tang_m')
        [EquinLT_backward, t] = AnEquin_all_forward_tang_m(L_prop_back, Equin_initial_prop_back', 0, 0, ...
            eps_max, x(end-14-14*(i-1)), x(end-13-14*(i-1)), 0, 0, 0, constants.mu, ...
            options.geopotential, options.third_body, constants.R_Earth, options.drag, 0);
    end
%     keyboard
    Equin_thrust = [Equin_thrust, EquinLT_backward];
    
    Equin_sorted{2*arcs-(i-1)*2} = EquinLT_backward(:,end:-1:1);
    Equin_sorted_ct(2*arcs-(i-1)*2) = 1;
    
    % Cost function: time spent using propulsion
    J = J - t(end);
%     if t(end)>0
%         keyboard
%     end
    % Time of flight update
    ToF = ToF - t(end);
    
    % Matching condition between final computed point of backward propagated propulsed
    % arc and required final point of backward propagated propulsed arc
    Ceq = [Ceq; EquinLT_backward(1:5,end) - Equin_final_prop_back(1:5)'];
    
    if options.gradient == 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gradient constraints
    index_GCeq =  [length(x) - 13 : length(x) - 12, length(x)-6 : length(x)] - (i-1)*14;
    eps_GCeq = eps_max;
    for i_GCeq = 1 : length(index_GCeq)
        x_var = x;
        
        if i_GCeq == length(x) - 6 -  (i-1)*14
            delta_temp = delta_a;
        else
            delta_temp = delta;
        end
        
        x_var(index_GCeq(i_GCeq)) = x_var(index_GCeq(i_GCeq)) + delta_temp;
        
%         [~,~,Ceq_var] = cost_constraint_function_OutOfPlane_2nd(x_var, parameters, options, constants);

        alpha_GCeq = x_var(end-13-14*(i-1));
        beta_GCeq  = x_var(end-12-14*(i-1));
        start_eq_GCeq = [x_var(end-5-14*(i-1):end-1-14*(i-1)) ...
            x_var(end-6-14*(i-1)) + x_var(end-14*(i-1))];
        start_eq_GCeq(1) = start_eq_GCeq(1) ;
        check_eq_GCeq = x_var(end-11-14*(i-1):end-6-14*(i-1))';
        check_eq_GCeq(1) = check_eq_GCeq(1) ;
        L_prop_GCeq = linspace(start_eq_GCeq(6), check_eq_GCeq(6), n);

        Ceq_var = jacobian_const_LT_OOP(eps_GCeq, alpha_GCeq, beta_GCeq, L_prop_GCeq, start_eq_GCeq, check_eq_GCeq(1:5), options, constants);
        Ceq_var(1) = Ceq_var(1);
        
        GCeq(temp,index_GCeq(i_GCeq)) = (Ceq_var(1) - Ceq(temp)) / delta_a;
        GCeq(temp+1:temp+4,index_GCeq(i_GCeq)) = (Ceq_var(2:end) - Ceq(temp+1:temp+4)) / delta;

        
    end
    
    GCeq(temp,   length(x) - 11 - (i-1)*14) = -1;
    GCeq(temp+1, length(x) - 10 - (i-1)*14) = -1;
    GCeq(temp+2, length(x) - 9 - (i-1)*14) = -1;
    GCeq(temp+3, length(x) - 8 - (i-1)*14) = -1;
    GCeq(temp+4, length(x) - 7 - (i-1)*14) = -1;
    
    temp = temp + 5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
    if i < arcs/2
        
        % -----------------------------------------------------------------
        % Initial conditions for next forward thrust arc
        % -----------------------------------------------------------------
        Equin_initial_prop_forw = x(17+14*(i-1): 22+14*(i-1));
        Equin_final_prop_forw   = [x(23+14*(i-1): 27+14*(i-1)) ...
            x(22+14*(i-1))+x(28+14*(i-1))];
        L_prop_forw = linspace(Equin_initial_prop_forw(6), Equin_final_prop_forw(6), n);
        
        
        % -----------------------------------------------------------------
        % Initial conditions for next backward thrust arc
        % -----------------------------------------------------------------
        Equin_initial_prop_back = [x(end-20-14*(i-1):end-16-14*(i-1))...
            x(end-21-14*(i-1)) + x(end-15-14*(i-1))];
        Equin_final_prop_back   = x(end-26-14*(i-1):end-21-14*(i-1));
        L_prop_back = linspace(Equin_initial_prop_back(6), Equin_final_prop_back(6), n);
        
        
        % ---------------------------------------------------------------------
        % Backward coast arc from final point of previous low thurst arc to initial point
        % of next propulsed arc
        % ---------------------------------------------------------------------
        
        % Backward coast arc
        L_coast_back = linspace(x(end-7-14*(i-1)), Equin_initial_prop_back(6), n);
        Equin_initial_coast_back = x(end-12-14*(i-1):end-7-14*(i-1));
        
%         [Equin_coast_backward, t] = AnEquin_all_forward_tang_m(L_coast_back, Equin_initial_coast_back, ...
%                                                                0, 0, ...
%                                                                0, 0, 0, ...
%                                                                0, 0, 0, ...
%                                                                constants.mu, options.geopotential, options.third_body,...
%                                                                constants.R_Earth, options.drag, 0);
    
%         [Equin_coast_backward, t] = AnEquin_forward_m(L_coast_back, Equin_initial_coast_back, 0, 0, 0, constants.mu);
        

        if strcmp(options.fun, 'AnEquin_forward_m')
            [Equin_coast_backward, t] = AnEquin_forward_m(L_coast_back, Equin_initial_coast_back, 0, 0, 0, constants.mu);
        elseif strcmp(options.fun, 'AnEquin_all_forward_tang_m')
            [Equin_coast_backward, t] = AnEquin_all_forward_tang_m(L_coast_back, Equin_initial_coast_back, 0, 0, ...
                0, 0, 0, 0, 0, 0, constants.mu, ...
                options.geopotential, options.third_body, constants.R_Earth, options.drag, 0);
        end

        Equin_coast = [Equin_coast, Equin_coast_backward];
        
        Equin_sorted{2*arcs-1-(i-1)*2} = Equin_coast_backward(:,end:-1:1);
        Equin_sorted_ct(2*arcs-1-(i-1)*2) = 0;
        
        % Matching condition between final point of backward coast arc and
        % initial point of next propulsed arc
        if size(Equin_coast_backward(1:5)' ) ~= size(x(end-20-14*(i-1):end-16-14*(i-1))') 
            Ceq = [Ceq; Equin_coast_backward(1:5) - x(end-20-14*(i-1):end-16-14*(i-1))'];
            
        else
            Ceq = [Ceq;Equin_coast_backward(1:5)' - x(end-20-14*(i-1):end-16-14*(i-1))'];
        end
        if options.gradient == 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gradient constraints
        index_GCeq =  [length(x) - 14, length(x)-11 : length(x)-6] - (i-1)*14;
        eps_GCeq = 0;
        for i_GCeq = 1 : length(index_GCeq)
            x_var = x;
            
            if i_GCeq == length(x) - 11 -  (i-1)*14
                delta_temp = delta_a;
            else
                delta_temp = delta;
            end
            
            x_var(index_GCeq(i_GCeq)) = x_var(index_GCeq(i_GCeq)) + delta_temp;
            
%             [~,~,Ceq_var] = cost_constraint_function_OutOfPlane_2nd(x_var, parameters, options, constants);
            
            alpha_GCeq = 0;
            beta_GCeq  = 0;
            start_eq_GCeq = x_var(end-11-14*(i-1):end-6-14*(i-1));
            start_eq_GCeq(1) =start_eq_GCeq(1) ;
            check_eq_GCeq = [x_var(end-19-14*(i-1):end-15-14*(i-1))';...
                x_var(end-20-14*(i-1))+x_var(end-14-14*(i-1))];
            check_eq_GCeq(1) =check_eq_GCeq(1) ;
            L_prop_GCeq = linspace(start_eq_GCeq(6), check_eq_GCeq(6), n);
            
            Ceq_var = jacobian_const_LT_OOP(eps_GCeq, alpha_GCeq, beta_GCeq, L_prop_GCeq, start_eq_GCeq, check_eq_GCeq(1:5), options, constants);
            Ceq_var(1) = Ceq_var(1) ;
            GCeq(temp,index_GCeq(i_GCeq)) = (Ceq_var(1) - Ceq(temp)) / delta_a;
            GCeq(temp+1:temp+4,index_GCeq(i_GCeq)) = (Ceq_var(2:end) - Ceq(temp+1:temp+4)) / delta;
        end
        
        GCeq(temp,   length(x) - 19 - (i-1)*14) = -1;
        GCeq(temp+1, length(x) - 18 - (i-1)*14) = -1;
        GCeq(temp+2, length(x) - 17 - (i-1)*14) = -1;
        GCeq(temp+3, length(x) - 16 - (i-1)*14) = -1;
        GCeq(temp+4, length(x) - 15 - (i-1)*14) = -1;
        
        
        temp = temp + 5;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        % Time of flight update
        ToF = ToF - t(end);
        
 
    end
    
    
    
end

if ~options.J_DeltaV && ~options.J_ToF
    J_out = 0;
elseif options.J_DeltaV 
    J_out = J;
elseif options.J_ToF
    J_out = ToF;
end

Ceq = [Ceq; (ToF - ToF_opt)*1e-3];


% Disequality non linear constraint
C = [x(end-7) + x(end-1) - arrival_eq(6); ...
    sum(x(14:14:end)) + departure_eq(6) - arrival_eq(6); -ToF];

% keyboard
% GCeq = sparse(GCeq);

%% Plot

if parameters.plot_flag
    
    % Thrusting and coasting arcs
    Cart_thrust = zeros(length(Equin_thrust),6);
    Cart_coast  = zeros(length(Equin_coast),6);
    for k = 1 : length(Equin_coast)
        Cart_coast(k,:) = eq2cart_m(Equin_coast(:,k),constants.mu);
    end
    for k = 1 : length(Equin_thrust)
        Cart_thrust(k,:) = eq2cart_m(Equin_thrust(:,k),constants.mu);
    end
    
    % Departure orbit
    departure_cart = eq2cart_m(departure_eq, constants.mu);
    
    % Eccentricity vector departure orbit
    departure_v = norm(departure_cart(4:6));
    departure_r = norm(departure_cart(1:3));
    departure_ecc_vector = 1 / constants.DU * ((departure_v^2 - constants.DU / departure_r) * departure_cart(1:3) - ...
                                            dot(departure_cart(1:3), departure_cart(4:6)) * departure_cart(4:6));
    
    % Arrival orbit
    arrival_cart = eq2cart_m(arrival_eq, constants.mu);

    % Eccentricity vector arrival orbit
    arrival_v = norm(arrival_cart(4:6));
    arrival_r = norm(arrival_cart(1:3));
    arrival_ecc_vector = (1 / constants.mu) * ((arrival_v^2 - constants.mu / arrival_r) * arrival_cart(1:3) - ...
                                            dot(arrival_cart(1:3), arrival_cart(4:6)) * arrival_cart(4:6));
                                        
    arrival_ecc_versor = arrival_ecc_vector ./ norm(arrival_ecc_vector);                                    
 
    % Period of departure and arrival trajectories
    T_departure = 2 * pi * sqrt(departure_eq(1)^3/constants.mu);
    T_arrival = 2 * pi * sqrt(arrival_eq(1)^3/constants.mu);
    
    figure
    hold on
    % Plot departure and arrival positions
    plot3(departure_cart(1), departure_cart(2), departure_cart(3),'bo','MarkerFaceColor','b','MarkerSize',8)
    plot3(arrival_cart(1), arrival_cart(2), arrival_cart(3),'ko','MarkerFaceColor','k','MarkerSize',8)
    
    % Plot departure and arrival trajectories
    plot_trajectory(departure_cart(1:3),departure_cart(4:6),...
        T_departure,constants.mu,'b',2);
    plot_trajectory(arrival_cart(1:3),arrival_cart(4:6),...
        T_arrival,constants.mu,'k',2);
    
    % Plot thrust and coast arcs
    index = 1;
    while index <= length(Equin_coast)
        plot3(Cart_coast(index : index + parameters.n - 2, 1),...
            Cart_coast(index : index + parameters.n - 2,2), ...
            Cart_coast(index : index + parameters.n - 2,3),'g-','LineWidth',2);
        index = index + parameters.n-1;
    end
    
    index = 1;
    while index <= length(Equin_thrust)
        plot3(Cart_thrust(index : index + parameters.n - 2, 1),...
            Cart_thrust(index : index + parameters.n - 2,2), ...
            Cart_thrust(index : index + parameters.n - 2,3),'r-','LineWidth',2);
        index = index + parameters.n-1;
    end
%     quiver3(0,0,0,1,0,0,'Color','b','LineWidth',2)
%     text(1,0,0,'Reference axis')
%     quiver3(0,0,0,0,1,0,'Color','b','LineWidth',2)
%     quiver3(0,0,0,arrival_ecc_versor(1),arrival_ecc_versor(2),arrival_ecc_versor(3),'Color','k','LineWidth',2)
%     text(arrival_ecc_versor(1),arrival_ecc_versor(2),arrival_ecc_versor(3),'Eccentricity vector final orbit')
%     quiver3(0,0,0,departure_ecc_vector(1),departure_ecc_vector(2),departure_ecc_vector(3),'Color','k','LineWidth',2)
    grid on
    xlabel('x [AU]')
    ylabel('y [AU]')
    zlabel('z [AU]')
    view(0,90)
    axis equal
    legend('Departure Point','Arrival Point','Initial orbit','Final orbit')
    
end







end




