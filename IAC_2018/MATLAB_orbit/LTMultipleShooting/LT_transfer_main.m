%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Low-thrust transfer between celestial bodies. Final arrival position is
% explicitly given.
% Inputs are: the keplerian or cartesian elements of the two objects
% (Departure and arrival) at the departure and arrival time, respectively

% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

clc
clear
addpath(genpath('../Propagator'))
addpath(genpath('../spaceart_toolbox'))

%% Constants

% Sun Gravitational Constant [AU^3/TU^2]
constants.mu = 1;

% Astronomical Unit [km]
constants.DU = 149597870.691;

% Time unit TU [days]
constants.TU = 58.13;

% Gravity acceleration [m/s2]
constants.g0 = 9.8;

% Seconds in a day
constants.sec_day = 86400;

% Adimensional Earth radius [DU]
constants.R_Earth = 6378.136 / constants.DU;


%% Low-thrust engine parameters

% Thrust [N]
engine.T = 0.3;

% Isp [s]
engine.Isp = 3000;

% Initial mass of the spacecraft[kg]
spacecraft.m = 2000;

% c [m/s]
engine.c = engine.Isp * constants.g0;

% m_dot [kg/s]
engine.m_dot = engine.T / engine.c;


%%  Input for transfer

spacecraft.m = 2000;

% Cartesian states at departure and arrival
state_dep = [2.3102068611528535, 1.4241695701267225, 0.1145892141646303, -0.2943515863842975, 0.5081468383730634, -0.030553944200837083];
state_arr = [-0.23169016737321624, 2.9095252345111184, 0.007231403980381036, -0.5844658829143383, 0.0109882200850904, 0.01894185269320052];

% Time of fligth [days]
ToF = 62319.3331828 - 62030.0;

% Adimensional time of flight
ToF = ToF / constants.TU;

%% Options for the optimisation

% Accelerations continuation method
options.k_acc = 1;


% Options for definition of P1 and P2 for initial guess: possible values
% are "departure" (same as departure point), "arrival" (same as arrival
% point), "linspace" (linear spaced from departure to arrival point) and
% "zeros" (All zero)
options.P10_P20 = 'departure';

% Impose coast arc before arrival?
options.final_coast = 0;

% Plot at the end?
options.plot_flag = 1;

% Impose constraint on time of flight
options.const_ToF = 1;

% Optimize deltaV?
options.J_DeltaV = 1;

% Nested functions for fmincon
options.nested =1;

% Give gradients to fmincon?
options.gradient = 0;

options.fun = 'AnEquin_forward_m';
% options.fun = 'AnEquin_all_forward_tang_m';

options.geopotential.J2 = 0;
options.geopotential.J3 = 0;
options.geopotential.J4 = 0;
options.geopotential.J5 = 0;

options.third_body.flag_Sun = 0;
options.third_body.flag_Moon = 0;

options.drag.CD = 0;

if options.gradient == 1
    options.options_fmincon_temp = optimset('Display','iter-detailed','MaxFunEvals',1e5,'LargeScale','off',...
        'Algorithm','interior-point','TolCon',5e-4 ,'TolX',1e-5, ...
        'GradConstr','off',...
        'GradObj','off',...
        'DerivativeCheck','off');
else
    options.options_fmincon_temp = optimset('Display','iter-detailed','MaxFunEvals',1e5,'LargeScale','off',...
        'Algorithm','interior-point','TolCon',5e-4 ,'TolX',1e-5, ...
        'GradConstr','off',...
        'GradObj','off');
end






%% Initialization and adimensionalization of variables

% Current acceleration for the transfer [m/s^2]
engine.acc = engine.T/ spacecraft.m;

% Initial keplerian elements
departure_kep = cart2kep(state_dep, constants.mu);

% Final keplerian elements
arrival_kep   = cart2kep(state_arr, constants.mu);

% Number of revolutions
n_rev =0;

options.n_rev = n_rev;

% Number of thrust arcs
options.transfer_arcs = 4 * (options.n_rev+1);

% Optimise transfer
[DeltaV, time, n_arcs_conv, xopt, Equin_sorted, Equin_sorted_ct, parameters] =...
    solve_NLP_LT(departure_kep, arrival_kep, ToF, engine, ....
    spacecraft, options, constants);

% Propellant mass
mf = spacecraft.m * exp(-DeltaV*1000/(engine.Isp * constants.g0));






