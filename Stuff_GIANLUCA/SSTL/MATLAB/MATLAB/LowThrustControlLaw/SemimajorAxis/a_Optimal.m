% =========================================================================
% Optimal thrusting angles strategy for change in the semimajor axis.
% Elevation equal to zero and azimuth in the direction of the flight path
% angle. 
% Both numerical and averaged analytical integration are realized in this
% file.
% Averaged analytical integration is performed with the analytical
% propagator of FABLE.
% User has to define the input (initial and final orbital elements,
% low-thrst acceleration).
% The script plot the variation of the orbital elements and return the
% DeltaV required (workspace)
% =========================================================================
% References: 
% A. E. Petropolous, "Simple Control Laws for low-thrust orbit
% transfers", AAS/AIAA Astrodynamics Specialist Conference, 2003
% A. Ruggiero et al. "Low-thrust maneuvers for the efficient correction of
% orbital elements", 32nd International Electric Propulsion Conference,
% 2011

% Marilena Di Carlo, 2015
% marilena.di-carlo@strath.ac.uk


clear
close all
addpath(genpath('..\..\spaceart_toolbox'))
addpath(genpath('..\..\Propagator'))
addpath('..\')


% Mode for the thrust control
mode = 'optimum_a';

% Targeted orbital element
element_target = 'a';

% =========================================================================
% Constants:
% =========================================================================
% Gravitational acceleration [km^3/s^2]
constants.mu_dim = 398600;
constants.mu_dim = 132712439935;

% Adimensional gravitational acceleration
constants.mu = 1;

% Time constants [s]
constants.TU = 806.78;
constants.TU = 58.13 * 86400;

% Length constant [km]
constants.DU = 6378.136;
constants.DU = 149597870.691;

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
% Input - user defined
% =========================================================================

% Spacecraft mass [kg]
m0 = 1000;

% Spacecraft thrust [N]
engine.T = 0.1;

% Acceleration [m/s^2]
engine.f = engine.T / m0;

% Acceleration [km/s^2]
engine.f = engine.f * 1e-3;

% Specific impulse [s]
engine.Isp = 1600;

% Initial orbital elements [km, rad]
initial_orbit.a      = 1 * constants.DU;
initial_orbit.e      = 0.1;
initial_orbit.i      = 30*pi/180;
initial_orbit.omega  = 60*pi/180;
initial_orbit.Omega  = 45*pi/180;
initial_orbit.E      = 180*pi/180;

% Final orbital elements [km, rad]
final_orbit.a =1.2*constants.DU;

if final_orbit.a > initial_orbit.a
elseif final_orbit.a < initial_orbit.a
    engine.f = -engine.f;
    engine.T = -engine.T;
elseif final_orbit.a == initial_orbit.a
    warning('Initial and final semiamjor axis are equal.');
    return
end

% Time for the integration - the integration stops when the desired element
% value is reached
time_int = 365 * 86400;

% =========================================================================
% Numerical Integration
% =========================================================================

% An event function is used to stop the integration when the targeted
% element has reached the desired value
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.a));
beta = 0;
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,engine.f,beta,constants, mode),[0 time_int], ...
                        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);

                    
Delta_V_numeric = Time(end) * abs(engine.f);                       
                    
% =========================================================================
% Averaged Analtyical Integration
% Here work with adimensional variables
% =========================================================================
[Delta_V_analytic, T_Analytic, Kep_Analytic] = semimajor_axis_adjustment(initial_orbit, final_orbit, time_int, engine, m0, element_target, constants);


% =========================================================================
% Analtyic equations
% =========================================================================
[K,E] = ellipke(4 * initial_orbit.e / (1 + initial_orbit.e)^2);
fe = 2 / (1 - initial_orbit.e) * E + 2 / (1 + initial_orbit.e) * K;
Delta_V_analytic_eq = 2 * pi * ...
    abs(sqrt(constants.mu_dim / initial_orbit.a) - sqrt(constants.mu_dim / final_orbit.a)) / ...
    (1 - initial_orbit.e^2) / fe;

a_analytic = ( 1 / initial_orbit.a + sign(engine.f) * Time.^2 * engine.f^2 * ...
    (1 - initial_orbit.e^2)^2 * fe^2 / (4 * pi^2 * constants.mu_dim) + ...
    -  engine.f * (1 - initial_orbit.e^2) * fe * Time / ...
    ( pi  * sqrt(constants.mu_dim * initial_orbit.a)) ).^(-1);

% =========================================================================
% Plots
% =========================================================================
% Numeric integration
figure('units','normalized','outerposition',[0 0 0.6 0.6])
subplot(2,2,1)
plot(Time/ 86400, Kep_Elements(:,1),'LineWidth',2,'Color','r')
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('Semimajor axis [km]')
subplot(2,2,2)
plot(Time/ 86400, Kep_Elements(:,2),'LineWidth',2,'Color','r')
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('Eccentricity')
subplot(2,3,4)
plot(Time/ 86400, Kep_Elements(:,3)*180/pi,'LineWidth',2,'Color','r')
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('Inclination [deg]')
subplot(2,3,6)
plot(Time/ 86400, Kep_Elements(:,5)*180/pi,'LineWidth',2,'Color','r')
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('\omega [deg]')  
subplot(2,3,5)
plot(Time/ 86400, Kep_Elements(:,4)*180/pi,'LineWidth',2,'Color','r')
grid on
xlim([0 Time(end)/86400])
xlabel('Time [days]')
ylabel('\Omega [deg]')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Numeric propagation', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

% Analytic integration
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
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Analytic propagation', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')



% Comparison numeric analytic
figure('units','normalized','outerposition',[0 0 0.6 0.6])
subplot(2,2,1)
hold on
plot(Time/ 86400, Kep_Elements(:,1),'LineWidth',2,'Color','r')
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,1) * constants.DU,'LineWidth',2,'Color','b')
legend('Numeric','Analytic')
grid on
xlim([0 T_Analytic(end) * constants.TU/86400])
xlabel('Time [days]')
ylabel('Semimajor axis [km]')
subplot(2,2,2)
hold on
plot(Time/ 86400, Kep_Elements(:,2),'LineWidth',2,'Color','r')
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,2),'LineWidth',2,'Color','b')
grid on
xlim([0 T_Analytic(end) * constants.TU/86400])
xlabel('Time [days]')
ylabel('Eccentricity')
subplot(2,3,4)
hold on
plot(Time/ 86400, Kep_Elements(:,3)*180/pi,'LineWidth',2,'Color','r')
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,3)*180/pi,'LineWidth',2,'Color','b')
grid on
xlim([0 T_Analytic(end)* constants.TU/86400])
xlabel('Time [days]')
ylabel('Inclination [deg]')
subplot(2,3,6)
hold on
plot(Time/ 86400, Kep_Elements(:,5)*180/pi,'LineWidth',2,'Color','r')
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,5)*180/pi,'LineWidth',2,'Color','b')
grid on
xlim([0 T_Analytic(end) * constants.TU/86400])
xlabel('Time [days]')
ylabel('\omega [deg]')  
subplot(2,3,5)
hold on
plot(Time/ 86400, Kep_Elements(:,4)*180/pi,'LineWidth',2,'Color','r')
plot(T_Analytic * constants.TU/ 86400, Kep_Analytic(:,4)*180/pi,'LineWidth',2,'Color','b')
grid on
xlim([0 T_Analytic(end) * constants.TU/86400])
xlabel('Time [days]')
ylabel('\Omega [deg]') 
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Comparison: numeric vs analytic', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
                     

% Analyic equations
figure
plot(Time/86400,a_analytic,'LineWidth',2)
xlabel('Time [days]')
ylabel('Semimajor axis [km]')
title('Analytic Equation')
grid on
