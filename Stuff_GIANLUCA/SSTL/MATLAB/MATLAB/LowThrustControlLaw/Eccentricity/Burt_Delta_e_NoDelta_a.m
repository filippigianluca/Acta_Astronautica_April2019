% =========================================================================
% Burt: change of eccentricity without change of semimajor axis
% Thrusting law to change the eccentricity without change of semimajor axis
% and argument of the perigee
% Numerical and analytical propagation (using FABLE)

% FABLE uses two thrust arc, with
% azimuth angle equal to 90 at the perigee and equal to 270 at the apogee.
% eta is zero. The averaged analytical integration is implemented
% but there is a problem for the semimajor axis. The semimajor axis value
% depends on the integral I11 that is not equal and opposite on the two
% thrust arcs because it is computed using P10 and P20, that are different
% at the beginning of the two thrust arcs. This causes a reduction of
% semimajor axis that should not be there. Investigate further. 
% The estimation of the deltaV however is correct.

% References: Burt "On space manoeuvres with continuous thrust"
% 
% Marilena Di Carlo, 2015
% marilena.di-carlo@strath.ac.uk
% =========================================================================

clear
close all
addpath(genpath('..\..\spaceart_toolbox'))
addpath('..\')
addpath(genpath('../../Propagator'))


mode = 'eccentricity_Delta_a=0';

element_target = 'eccentricity';

% =========================================================================
% Constants:
% =========================================================================
% Gravitational acceleration [km^3/s^2]
constants.mu_dim = 398600;

% Adimensional gravitational acceleration
constants.mu = 1;

% Time constants [s]
constants.TU = 806.78;

% Length constant [km]
constants.DU = 6378.136;

% Gravity acceleration [m/s^2]
constants.g0_m_s = 9.8;

% Gravity acceleration [km/s^2]
constants.g0_km_s = constants.g0_m_s * 1e-3;

% Gravity acceleration [DU/TU^2]
constants.g0 = constants.g0_m_s * 1e-3 * constants.TU^2 / constants.DU;

% No J2 considered
constants.J2 = 0;

% Earth radius [km]
constants.R_Earth_dim = 6378;
constants.R_Earth = 1;


% =========================================================================
% Input - defined by the user
% =========================================================================

% Spacecraft mass [kg]
m0 = 2000;

% Spacecraft thrust [N]
engine.T = 0.3;

% Acceleration [m/s^2]
engine.f = engine.T / m0;

% Acceleration [km/s^2]
engine.f = engine.f * 1e-3;

% Specific impulse [s]
engine.Isp = 3000;

% Specific impulse [TU]
% Isp = Isp / constants.TU;

% Initial orbit elements [km, rad]
initial_orbit.a     = 7000;
initial_orbit.e     = 0.005;
initial_orbit.i     = 30*pi/180;
initial_orbit.omega = 00*pi/180;
initial_orbit.Omega = 0*pi/180;
initial_orbit.E = 10*pi/180;

% Final orbit [km, rad]
final_orbit.e = 0.05;

if final_orbit.e < initial_orbit.e
    engine.f = -engine.f;
    k_p = -1;
    k_a = 1;
elseif final_orbit.e > initial_orbit.e
    k_p = 1;
    k_a = -1;
elseif final_orbit.e == initial_orbit.e
    warning('Initial and final eccentricity are equal.');
end

time_int =365*86400;

% =========================================================================
% Numerical Integration
% =========================================================================

% An event function is used to stop the integration when the targeted
% element has reached the desired value

beta = 0;
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.e));
tic
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 time_int], ...
                        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);
toc
     
Delta_V_numeric = Time(end) * abs(engine.f);


% =========================================================================
% Analtyic equations
% =========================================================================
Delta_V_equation = pi/4 * sqrt(constants.mu_dim / initial_orbit.a) * abs(asin(initial_orbit.e) - asin(final_orbit.e));    

e_equation = sin(4/pi * sqrt(initial_orbit.a/constants.mu_dim) * engine.f * Time + asin(initial_orbit.e));

% =========================================================================
% Averaged Analtyical Integration
% Here work with DU and TU
% =========================================================================

[Delta_V_analytic, T_Analytic, Kep_Analytic] = Burt_a0_e(initial_orbit, final_orbit, time_int, engine, m0, element_target, constants);


% =========================================================================
% Plots
% =========================================================================
figure('units','normalized','outerposition',[0 0 0.6 0.6])
subplot(2,2,1)
plot(Time/ 86400, Kep_Elements(:,1),'LineWidth',2)
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('Semimajor axis [km]')
subplot(2,2,2)
plot(Time/ 86400, Kep_Elements(:,2),'LineWidth',2)
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('Eccentricity')
subplot(2,3,4)
plot(Time/ 86400, Kep_Elements(:,3)*180/pi,'LineWidth',2)
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('Inclination [deg]')
subplot(2,3,6)
plot(Time/ 86400, Kep_Elements(:,5)*180/pi,'LineWidth',2)
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('\omega [deg]')  
subplot(2,3,5)
plot(Time/ 86400, Kep_Elements(:,4)*180/pi,'LineWidth',2)
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('\Omega [deg]')


figure('units','normalized','outerposition',[0 0 0.6 0.6])
subplot(2,2,1)
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,1) * constants.DU,'LineWidth',2)
grid on
xlim([0 T_Analytic(end) * constants.TU/86400])
xlabel('Time [days]')
ylabel('Semimajor axis [km]')
subplot(2,2,2)
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,2),'LineWidth',2)
grid on
xlim([0 T_Analytic(end) * constants.TU/86400])
xlabel('Time [days]')
ylabel('Eccentricity')
subplot(2,3,4)
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,3)*180/pi,'LineWidth',2)
grid on
xlim([0 T_Analytic(end)* constants.TU/86400])
xlabel('Time [days]')
ylabel('Inclination [deg]')
subplot(2,3,6)
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,5)*180/pi,'LineWidth',2)
grid on
xlim([0 T_Analytic(end) * constants.TU/86400])
xlabel('Time [days]')
ylabel('\omega [deg]')  
subplot(2,3,5)
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,4)*180/pi,'LineWidth',2)
grid on
xlim([0 T_Analytic(end) * constants.TU/86400])
xlabel('Time [days]')
ylabel('\Omega [deg]')  
                     


% Comparison
figure('units','normalized','outerposition',[0 0 0.6 0.6])
subplot(2,2,1)
hold on
plot(Time/ 86400, Kep_Elements(:,1),'LineWidth',2,'Color','r')
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,1) * constants.DU,'LineWidth',2,'Color','b')
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('Semimajor axis [km]')
legend('Numerical','Analytic')
subplot(2,2,2)
hold on
plot(Time/ 86400, Kep_Elements(:,2),'LineWidth',2,'Color','r')
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,2) ,'LineWidth',2,'Color','b')
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('Eccentricity')
subplot(2,3,4)
hold on
plot(Time/ 86400, Kep_Elements(:,3)*180/pi,'LineWidth',2,'Color','r')
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,3) *180/pi,'LineWidth',2,'Color','b')
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('Inclination [deg]')
subplot(2,3,6)
hold on
plot(Time/ 86400, Kep_Elements(:,5)*180/pi,'LineWidth',2,'Color','r')
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,5) * 180/pi,'LineWidth',2,'Color','b')
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('\omega [deg]')  
subplot(2,3,5)
hold on
plot(Time/ 86400, Kep_Elements(:,4)*180/pi,'LineWidth',2,'Color','r')
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,4)*180/pi,'LineWidth',2,'Color','b')
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('\Omega [deg]')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Comparison', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')


figure
plot(Time/86400, e_equation, 'LineWidth',2)
xlabel('Time [days]')
ylabel('e')
title('From analytic equation')
grid on
