%==========================================================================
% Variation of semimajor axis, inclination and right ascension with J2 in a
% given time of flight using low-thrust - Third strategy
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

% Total ToF [s]
ToF_total =600*86400;


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

tic
[DeltaV, alpha1, ToF1, alpha2, ToF2] = aRAAN_i_2arcs(initial_orbit, final_orbit, ToF_total, ...
                                    f, constants);
toc

%% Plot

% -------------------------------------------------------------------------
% First phase
% -------------------------------------------------------------------------
% Vector of time for the first phase
t1 = linspace(0, ToF1, 1000);

% Velocity variation during the first phase
Vt1 = sqrt(constants.mu/a_initial) - 2 * f * alpha1 * t1 / pi;

% Variation of semimajor axis during the first phase
at1 = constants.mu ./ Vt1.^2;

% Variation of right ascension during the first phase
RAAN_t1 = mod(RAAN_initial + ...
    3/32 * pi * constants.mu * constants.J2 * constants.R_Earth^2 * cos(i_initial) / (f * alpha1) * (1./at1.^4 - 1/a_initial^4), 2*pi);

% -------------------------------------------------------------------------
% Second phase
% -------------------------------------------------------------------------
% Time vector
t2 = linspace(0, ToF2, 1000);

% Inclination
it2 = i_initial + 2 * f * a_final * sin(alpha2) / (pi * sqrt(constants.mu * a_final)) * t2;

% Right ascension
RAAN_t2 = mod(RAAN_t1(end) + ...
    3/4 * pi * constants.mu * constants.J2 * constants.R_Earth^2  / (a_final^4 * f* sin(alpha2)) * (sin(i_initial) - sin(it2)), 2*pi);


% Plot semimajor axis
ls = 16;
figure(1)
hold on
plot(t1/86400,at1,'LineWidth',2)
plot(t1(end)/86400 + t2/86400, at1(end)*ones(1,length(t2)),'LineWidth',2)
grid on
xlabel('Time  [days]','FontSize',ls)
ylabel('Semimajor axis [km]','FontSize',ls)
set(gca,'fontsize',ls)

% Plot right ascension
figure(2)
hold on
plot(t1/86400,RAAN_t1*180/pi,'LineWidth',2)
plot(t1(end)/86400 + t2/86400,RAAN_t2*180/pi,'LineWidth',2)
grid on
xlabel('Time  [days]','FontSize',ls)
ylabel('RAAN [deg]','FontSize',ls)
set(gca,'fontsize',ls)

% Plot inclination
figure(3)
hold on
plot(t1/86400,i_initial * 180/ pi * ones(1,length(t1)),'LineWidth',2)
plot(t1(end)/86400 + t2/86400,it2*180/pi,'LineWidth',2)
grid on
xlabel('Time [days]','FontSize',ls)
ylabel('Inclination [deg]','FontSize',ls)
set(gca,'fontsize',ls)
ylim([50 56])




