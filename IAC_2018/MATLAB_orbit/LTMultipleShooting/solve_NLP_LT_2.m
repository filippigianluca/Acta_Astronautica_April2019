function [DeltaV, time, number_arcs_convergence, ...
    xopt, Equin_sorted, Equin_sorted_ct, parameters] = solve_NLP_LT_2 ...
    (departure_kep, arrival_kep, engine, spacecraft, options, constants)


% Function for solution of planar low thrust transfer problem
% Input: departure_kep -> keplerian elements of the departure point at time
% t0
%        arrival_kep   -> keplerian elements of the arrival point at time
%        t0
%        engine        -> structure with information about the low-thrust
%                         engine (thrust, acceleration, specific impulse)
%        spacecraft    -> structure with information about the spacecrat
%                         (could be removed)
%        options       -> options for the transfer (number of transfer
%                         arcs, continuation method...)
%        constants     -> structure with constants

% Output: DeltaV       -> Delta V for the transfer
%         xopt         ->
%
% Author: Marilena Di Carlo, 2017
% email: marilena.di-carlo@strath.ac.uk



%%

% Equinoctial elements of the DEPARTURE point [AU, rad]
departure_eq = kep2eq(departure_kep);


% new implementation: number of arcs defined by the user
arcs = options.transfer_arcs;


eps_max = engine.acc;

% Low-thrust acceleration [AU/TU^2]
eps_max = eps_max * 10^(-3) * (((constants.TU*constants.sec_day)^2) / constants.DU);


% =====================================================================
% Generation of nodes along the orbit for the starting point of thrust and coast arcs
% =====================================================================

% Variation of true longitude from departure to arrival point
dL = linspace(departure_eq(6), departure_eq(6)+2*pi, 2*arcs+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
equin_nodes(1,:) = linspace(departure_eq(1), arrival_kep(1), 2*arcs);
equin_nodes(2,:) = departure_eq(2) * ones(1,2*arcs);
equin_nodes(3,:) = departure_eq(3) * ones(1,2*arcs);
equin_nodes(4,:) = departure_eq(4) * ones(1,2*arcs);
equin_nodes(5,:) = departure_eq(5) * ones(1,2*arcs);
for i = 1 : size(equin_nodes,2)
    equin_nodes(6,i) = dL(i+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% =======================================================
% Initial guess vector
% ========================================================

% The vector x is defined by 10 elements per arc.
% x(1)     - azimuth angle for thrust vector
% x(2)      - elevation angle
% x(3:6)   - orbital elements of the initial point of the 1st thrust arc
%           (a, h, k, L)
% x(7:10)  - orbital elements of the intial poinf of the coast arc (a, h, k, L)
% ... and so on

x00 = zeros(arcs, 14);

% Azimuth angle
if arrival_kep(1) >= departure_eq(1)
    x00(:,1) =90*(pi/180);
else
    x00(:,1) = -90*(pi/180);
end
x00(:,2) = 0*(pi/180);

for i = 1 :  arcs
    
    x00(i,3:8)  = [equin_nodes(1:5,1+(i-1)*2)' equin_nodes(6,1+(i-1)*2)];
    x00(i,9:14)  = [equin_nodes(1:5,2*i)' equin_nodes(6,2*i)-equin_nodes(6,1+(i-1)*2)];
    
end


x0 = reshape(x00', 1, 14*arcs);

% Initial guess for time of flight
if options.flag_IG == 0 
    x0 = [x0, options.IG_ToF / constants.TU];
elseif options.flag_IG
    x0 = options.IG;
end

% =========================================================
% Lower and upper boundary for the states
% =========================================================
min_semimajor = min(departure_eq(1), arrival_kep(1));
max_semimajor = max(departure_eq(1), arrival_kep(1));


LB = reshape([-180* (pi/180);  -90*pi/180;   min_semimajor - 0.5 * min_semimajor;  -1; -1;   -Inf; -Inf; departure_eq(6);  ...
    min_semimajor - 0.5 * min_semimajor;  -1; -1;   -Inf; -Inf; 0] * ones(1,arcs) , 1, arcs+13*arcs);
UB = reshape([ 180* (pi/180);   90*pi/180;   max_semimajor + 0.5 * max_semimajor;   1;  1;   Inf; Inf; Inf; ...
    max_semimajor + 0.5 * max_semimajor;   1;  1;   Inf; Inf; Inf]   * ones(1,arcs) , 1, arcs+13*arcs);


% Lower and upper boundaries for the time of fligth
LB = [LB, options.LB_ToF/constants.TU];
UB = [UB, options.UB_ToF/constants.TU];


% =========================================================
% Linear constraints
% =========================================================
%
% INEQUALITIES CONSTRAINTS ON TRUE LONGITUDE VALUES
% Inequalities conditions:
% departure_eq(6) < L1 < L2 < L3 and so on until Ln < arrival_eq(6)
% Inequalities constraints are expressed through the relationship Ax<b

A = zeros(2*arcs-2, 14*arcs+1);
b = zeros(2*arcs-2,1);
for i = 1 : arcs-1
    A(i ,8 + 14 * (i-1)) = 1;
    A(i, 8+14+ 14 * (i-1)) = -1;
end
for i = 1 : arcs-1
    A(arcs-1+i,8 + 14 * (i-1)) = 1;
    A(arcs-1+i,14 + 14 * (i-1)) = 1;
    A(arcs-1+i,8 +14+ 14 * (i-1)) = -1;
end
% A(end-1,8 + 14 * (i)) = 1;
% A(end-1,14 + 14 * (i)) = 1;
% b(end-1,1) = arrival_eq(6);
% A(end,14:14:end) = 1;
% b(end) = 2*pi*options.n_rev;
% b(end) = arrival_eq(6) - departure_eq(6);




% =========================================================
% Optimization
% =========================================================

parameters.arcs         = arcs;
parameters.eps_max      = eps_max;
parameters.departure_eq = departure_eq;
parameters.arrival_kep_t0 = arrival_kep;
parameters.plot_flag    = 0;
parameters.n            = 2;
parameters.m            = spacecraft.m;


[~,~,Ceq_x0,~,~,~,~] = cost_constraint_function_LT_2(x0, parameters, options, constants);
x0_C_viol = max(abs(Ceq_x0));

if x0_C_viol == 0
    DeltaV = 0;
    time = 0;
    number_arcs_convergence = 0;
    xopt = [];
    Equin_sorted = [];
    Equin_sorted_ct=[];
else
    
    options.options_fmincon = optimset('Display',options.options_fmincon_temp.Display,...
        'MaxFunEvals',options.options_fmincon_temp.MaxFunEvals,...
        'LargeScale',options.options_fmincon_temp.LargeScale,...
        'Algorithm',options.options_fmincon_temp.Algorithm,...
        'TolCon', options.options_fmincon_temp.TolCon /x0_C_viol,...
        'TolX',options.options_fmincon_temp.TolX, ...
        'GradConstr',options.options_fmincon_temp.GradConstr,...
        'GradObj',options.options_fmincon_temp.GradObj,...
        'DerivativeCheck','off');

        tic
        [xopt,fopt,iout,output] = runobjconstr_LT_2(x0, parameters, options, ...
            constants, A, b, LB, UB);
        toc
    parameters.plot_flag = options.plot_flag;
    parameters.n         = 100;
    
    [~,C_opt,Ceq_opt2,~,ToF,Equin_sorted,Equin_sorted_ct,J_DV] = cost_constraint_function_LT_2(xopt, parameters, options, constants);
    max(abs(Ceq_opt2))
    time = ToF * constants.TU;
    
    
    
    if (iout == 1  ||  iout == 2)
        
        % deltaV in km/s (AU is in km)
        DeltaV = J_DV * eps_max / ((constants.sec_day*constants.TU)/constants.DU);
        
        % Time spent using low thrust propulsion
        LowThrustTime = J_DV * constants.TU;
        
        number_arcs_convergence = options.transfer_arcs;
        
       
    else
        DeltaV = NaN;
        number_arcs_convergence = NaN;
        LowThrustTime = NaN;
    end
  
end
end

