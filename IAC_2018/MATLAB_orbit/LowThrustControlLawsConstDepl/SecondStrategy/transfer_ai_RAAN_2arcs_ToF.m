%==========================================================================
% Variation of semimajor axis, inclination and right ascension with J2 in a
% given time of flight using low-thrust - Second strategy
% =========================================================================
% Reference: M. Di Carlo, A. Ricciardi, M. Vasile, "Optimised Constellation
% Deployment using Low-Thrust Propulsion", SPACE 2016
% Marilena Di Carlo, 2016
% marilena.di-carlo@strath.ac.uk


clear 
close all
clc

addpath(genpath('../'))

%% Constants

% Gravitational constant of the Earth [km^3/s^2]
constants.mu = 398600;

% J2
constants.J2 = 1.0826e-3;

% Earth radius [km]
constants.R_Earth = 6378.136;



%% Input

% Initial semimajor axis [km]
a_initial = 10000;
% Final semimajor axis [km[
a_final   = 24200;

% Initial inclination
i_initial = 51* pi/180;
% Final inclination
i_final = 56 * pi/180;

% Initial RAAN
RAAN_initial = 360 * pi/180;
% Final RAAN
RAAN_final   = 150 * pi/180;

% ?
f_direction = 1;

% Total ToF [s]
ToF_total = 600*86400;


% Low thrust acceleration [km/s^2]
f = 1.5e-7;



%% Initialise variables

% Initial orbit
initial_orbit.a    = a_initial;
initial_orbit.incl = i_initial;
initial_orbit.RAAN = RAAN_initial;

% Final orbit
final_orbit.a    = a_final;
final_orbit.incl = i_final;
final_orbit.RAAN = RAAN_final;


%% Compute transfer

% Since number of equations is lower than number of unknowns, the transfer
% is attempted with different ToF for the first phase. 
% all_solutions collect all the results.
% best_solution gives detail about the minimum deltaV result

[best_solution, all_solutions] = ai_RAAN_2arcs(initial_orbit, final_orbit, ...
                                    ToF_total, f_direction, ...
                                    f, constants);
       
% -------------------------------------------------------------------------
% Minimum DeltaV solution
% -------------------------------------------------------------------------
% DeltaV minimum DeltaV solution
DeltaV     = best_solution.DeltaV;
% Length thrust arc first phase
alpha1_min = best_solution.alpha1_min;
% Elevation angle first phase
beta       = best_solution.beta;
% Lenght thrust arc second phase
alpha2_min = best_solution.alpha2_min;
% Time of flight first and second phase
ToF1st     = best_solution.T1st_min;
ToF2nd     = best_solution.T2nd_min;


% -------------------------------------------------------------------------
% ALL solutions
% -------------------------------------------------------------------------
% All DeltaV
DeltaV_TOT = all_solutions.DeltaV_TOT;
% All length thrust arc first phase
alpha1     = all_solutions.alpha1;
% All lenght thrust arc second phase
alpha2     = all_solutions.alpha2;




%% Plot

% Time interval first phase
t1 = linspace(0, ToF1st, 1000);

% Velocity variation during the first phase
% Using alpha1_min as value for the length of the thrust arcs and beta as
% value for the elevation angle
V_t = sqrt(constants.mu/a_initial) - 2 * f * cos(beta) * alpha1_min * t1 / pi;

% Semimajor axis variation during the first phase
a_t = constants.mu./V_t.^2;

% Inclination variation during the first phsase
i_t = i_initial + tan(beta) * sin(alpha1_min) / alpha1_min * log(a_t/a_initial)/2;

% Final right ascension at the end of the first phase
plot_flag = 0;
[~, ~, RAAN_ToF] = a_i_circular_2arcs(initial_orbit, final_orbit, ...
    alpha1_min, f, constants, plot_flag);

% Rate of change of RAAN during the second phase due to J2
RAAN_dot = 3/2* sqrt(constants.mu) * constants.J2 * constants.R_Earth^2 * cos(i_final) * ...
    a_final^(-7/2);

% Time interval of the second phase
t2 = linspace(0, ToF2nd, 1000);

% Variation of right ascension during the second phase: due to J2 and low
% thrust engine
for i = 1 : length(t2)
    RAAN_t(i) = mod(RAAN_ToF + ...
        f_direction * 2 * f * sin(alpha2_min)/ (pi * sin(i_final)) * sqrt(a_final/constants.mu)  * t2(i)+ ...
        - RAAN_dot * t2(i), 2*pi);
end

% Variation of right ascension during the first phase
k2 = - 4 * log(a_final/a_initial) / (i_final - i_initial);

k3 = exp(k2 * i_t) .* (k2 * cos(i_t) + sin(i_t)) - ...
    exp(k2 * i_initial) * (k2 * cos(i_initial) + sin(i_initial));

k1 = -3/4 * constants.mu * pi * constants.J2 * constants.R_Earth^2 / (a_initial^4 * f * sin(beta) * sin(alpha1_min)) * ...
    exp(4 * log(a_final/a_initial) * i_initial / (i_final - i_initial));


RAAN_0 = mod(RAAN_initial + k1 / (1+k2^2) * k3, 2*pi);


% -------------------------------------------------------------------------
% Plot semimajor axis
% -------------------------------------------------------------------------
ls = 16;
figure(1)
hold on
% variation of semiajor axis during the first phase
plot(t1/86400, a_t,'LineWidth',2)
% Constant semimajor axis during the second phase
plot(t1(end)/86400 + t2/86400, a_final * ones(1,length(t2)),'LineWidth',2)
grid on
xlabel('Time [days]','FontSize',ls)
ylabel('Semimajor axis [km]','FontSize',ls)
set(gca,'fontsize',ls)


% -------------------------------------------------------------------------
% Plot inclination
% -------------------------------------------------------------------------
figure(2)
hold on
% Variation of inclination during the first phase
plot(t1/86400, i_t*180/pi,'LineWidth',2)
% Constant inclination during the second phase
plot(t1(end)/86400 + t2/86400, i_t(end)*180/pi* ones(1,length(t2)),'LineWidth',2)
grid on
xlabel('Time [days]','FontSize',ls)
ylabel('Inclination [deg]','FontSize',ls)
set(gca,'fontsize',ls)
ylim([50 56])

% -------------------------------------------------------------------------
% Plot right ascension
% -------------------------------------------------------------------------
figure(3)
hold on
% Variation of right ascension during the first phase
plot(t1/86400, RAAN_0*180/pi,'LineWidth',2)
% Variation of right ascnesion dring the second phase
plot(t1(end)/86400 + t2/86400, RAAN_t*180/pi,'LineWidth',2)
hold on
grid on
xlabel('Time [days]','FontSize',ls)
ylabel('\Omega [deg]','FontSize',ls)
set(gca,'fontsize',ls)

            
% ----------------------------------------------------------
figure
plot(all_solutions.ToF_1st,...
    all_solutions.DeltaV_TOT)
plot(all_solutions.ToF_1st/86400,all_solutions.DeltaV_TOT,'o')
grid on
xlabel('ToF_1 [days]','FontSize',ls)
ylabel('\Delta V [km/s]','FontSize',ls)
set(gca,'FontSize',ls)