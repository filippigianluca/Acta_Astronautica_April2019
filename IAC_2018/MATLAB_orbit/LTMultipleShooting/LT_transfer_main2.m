%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Low-thrust transfer between celestial bodies. Final arrival position is
% not given, but computed based on transfer time.
% Direct optimisation with multiple shooting and analytic propagation.
% Problems that can be solved:
% - feasibility (go from initial to final body in any time of flight and
% with any deltaV)
% - minimise time of fligth (go from initial to final body with minimum
% time of flight)
% - minimise DeltaV
% At the moment it does not handle transfer longer than one orbital period.
% Fix this
% Inputs are: the keplerian or cartesian elements of the two objects
% (Departure and arrival) at the same epoch (initial time)

% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk


addpath(genpath('../Propagator'))
addpath(genpath('../spaceart_toolbox'))




%% Constants

% Sun Gravitational Constant [AU^3/TU^2]
constants.mu = 1.0;

% Astronomical Unit [km]
constants.DU = 149597870.691;

% Gravity acceleration [m/s2]
constants.g0 = 9.80665;

% Seconds in a day
constants.sec_day = 86400;

% Time unit TU [days]
constants.TU = 5022642.89091/constants.sec_day;

% Adimensional Earth radius [DU]
constants.R_Earth = 6378.136 / constants.DU;



%% INPUT

% Spacecraft mass
spacecraft.m = 1000;

% Keplerian elements of departure and arrival body at the same initial time
% (initial time is time of departure)
kep_dep = [2.3614601	0.0886122	7.1404900*pi/180	151.2172200*pi/180	103.8512900*pi/180	326.5320247*pi/180];
kep_arr0 = [2.5	0.0886122	7.1404900*pi/180	151.2172200*pi/180	103.8512900*pi/180	330*pi/180];

% Cartesian elements of bodies at departure
state_dep = kep2cart(kep_dep, constants.mu);
state_arr0 =  kep2cart(kep_arr0, constants.mu);

%% Low-thrust engine parameters

% Thrust [N]
engine.T = 0.3;

% Isp [s]
engine.Isp = 3000.0;

% c [m/s]
engine.c = engine.Isp * constants.g0;

% m_dot [kg/s]
engine.m_dot = engine.T / engine.c;



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

% ----------------- What do you want to minimise? -------------------------
% Optimize deltaV?
options.J_DeltaV = 0;

% Optimise time of flight?
options.J_ToF = 1;
% -------------------------------------------------------------------------

% Give gradients to fmincon?
options.gradient = 0;

% Initial guess?
options.flag_IG  = 0;
% options.IG = initial_guess;

% Lower and upper values for time of flight [days] - If you want transfer
% in a given time of flight, set LB equal to UB
options.LB_ToF = 100;
options.UB_ToF = 1000;

% Initial guess time of flight [days]
options.IG_ToF = 200;

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
    options.options_fmincon_temp = optimset('Display','off','MaxFunEvals',1e5,'LargeScale','off',...
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
arrival_kep   = cart2kep(state_arr0, constants.mu);

% Number of revolutions
n_rev =0;

options.n_rev = n_rev;

% Number of thrust arcs
options.transfer_arcs = 4 * (options.n_rev+1);

%%  Optimise transfer - minimisation of DeltaV for transfer in given ToF


[DeltaV, time, n_arcs_conv, xopt, Equin_sorted, Equin_sorted_ct, parameters] =...
    solve_NLP_LT_2(departure_kep, arrival_kep,  engine, ....
    spacecraft, options, constants);

% Propellant mass
mf_out = spacecraft.m * exp(-DeltaV*1000/(engine.Isp * constants.g0));

tof_out = time;
converged = n_arcs_conv == options.transfer_arcs;

