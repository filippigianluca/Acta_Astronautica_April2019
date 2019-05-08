% =========================================================================
% Burt: change of omega without change of semimajor axis and eccentricity
% Numerical and analytical propagation (using FABLE).
% Second case. 
% FABLE computes an integration with two thrust arcs, with the first thrust
% arc centered at 90 deg. The acceleration is only tangential and changes
% sign on the two arcs.

% References: Burt "On space manoeuvres with continuous thrust"
% 
% Marilena Di Carlo, 2016
% marilena.di-carlo@strath.ac.uk
% =========================================================================

clear
close all
addpath(genpath('..\..\spaceart_toolbox'))
addpath('..\')
addpath(genpath('..\..\Propagator'))


mode = 'omega_Delta_a_Delta_e=0_2';

element_target = 'omega';

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
m0 = 3000;

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
initial_orbit.e     = 0.1;
initial_orbit.i     = 30*pi/180;
initial_orbit.omega = 0*pi/180;
initial_orbit.Omega = 30*pi/180;
initial_orbit.E = 10*pi/180;

% Final orbit [km, rad]
final_orbit.omega = 20 * pi/180;

if final_orbit.omega < initial_orbit.omega
    engine.f = -engine.f;
%     k_p = -1;
%     k_a = 1;
elseif final_orbit.omega > initial_orbit.omega
%     k_p = 1;
%     k_a = -1;
elseif final_orbit.omega == initial_orbit.omega
    warning('Initial and final omega are equal.');
end

time_int =365*86400;

% =========================================================================
% Numerical Integration
% =========================================================================

% An event function is used to stop the integration when the targeted
% element has reached the desired value

beta = 0;
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.omega));
tic
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 time_int], ...
                        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);
toc
     
Delta_V_numeric = Time(end) * abs(engine.f);



% =========================================================================
% Averaged Analtyical Integration
% Here work with DU and TU
% =========================================================================
Delta_V_equation = abs(final_orbit.omega - initial_orbit.omega) / ...
    (2 / pi * sqrt(initial_orbit.a / constants.mu_dim) * (2 - initial_orbit.e^2) / initial_orbit.e);

omega_analytic = initial_orbit.omega + 2/pi * sqrt(initial_orbit.a / constants.mu_dim) * ...
    (2 - initial_orbit.e^2) / initial_orbit.e * engine.T/m0 * 1e-3 * Time;

[Delta_V_analytic, T_Analytic, Kep_Analytic] = Burt_a0_e0_omega2(initial_orbit, final_orbit, time_int, engine, m0, element_target, constants);


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


% =========================================================================
% Analytic equations
% =========================================================================
figure
plot(Time/86400, omega_analytic*180/pi,'LineWidth',2)
xlabel('Time [days]')
ylabel('\omega [deg]')  
grid on
title('Analytic equation')
