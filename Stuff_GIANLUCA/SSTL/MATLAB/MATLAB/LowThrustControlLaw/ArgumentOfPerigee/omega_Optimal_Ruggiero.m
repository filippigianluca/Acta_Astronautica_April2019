% =========================================================================
% Optimal thrusting angles strategy for change in the argument of perigee.
% Reference: Ruggiero et al., "Low-thrust manuevers for the efficient
% correction of orbital elements"
% Numerical integration of the orbital elements
% =========================================================================
% Marilena Di Carlo, 2015
% marilena.di-carlo@strath.ac.uk



clear
close all
addpath(genpath('spaceart_toolbox'))
addpath('../')

% Mode for the thrust control
mode = 'optimum_omega_Ruggiero';


% Targeted orbital element
element_target = 'omega';

% =========================================================================
% Constants:
% =========================================================================
% Gravitational acceleration [km^3/s^2]
constants.mu_dim = 398600;

constants.mu = 1;

% Time constants [s]
constants.TU = 806.78;

% Length constant [km]
constants.DU = 6378.136;

% Gravity acceleration [m/s^2]
constants.g0 = 9.8;

% Gravity acceleration [DU/TU^2]
constants.g0_DUTU = constants.g0 * 1e-3 * constants.TU^2 / constants.DU;


% =========================================================================
% Input - user defined
% =========================================================================

% Spacecraft mass [kg]
m0 = 3000;

% Spacecraft thrust [N]
T = 0.3;

% Acceleration [m/s^2]
f = T / m0;

% Acceleration [km/s^2]
f = f * 1e-3;

% Specific impulse [s]
Isp = 3000;

% Specific impulse [TU]
% Isp = Isp / constants.TU;

% Initial orbit elements [km, rad]
initial_orbit.a     = 7000;
initial_orbit.e     = 0.1;
initial_orbit.i     = 30*pi/180;
initial_orbit.omega = 0*pi/180;
initial_orbit.Omega  = 30*pi/180;
initial_orbit.E = 80*pi/180;

% Final orbit [km, rad]
final_orbit.omega = 40*pi/180;

% Time for the integration - the integration stops when the desired element
% value is reached
time_int = 300 * 86400;

% =========================================================================
% Numerical Integration
% =========================================================================

% An event function is used to stop the integration when the targeted
% element has reached the desired value

if final_orbit.omega < initial_orbit.omega
    f = -f;
end
options_ODE = odeset('AbsTol',1e-8,'RelTol',1e-8,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.omega));
tic
beta = 0;
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 time_int], ...
                        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);
                  
Delta_V_numeric = Time(end) * abs(f);


% =========================================================================
% Analytic equations
% =========================================================================
omega_analytic = initial_orbit.omega + sqrt(initial_orbit.a / constants.mu_dim) * ...
    sqrt(1 - initial_orbit.e^2) * (3 - 2 * initial_orbit.e^2) * f / ...
    (2 * initial_orbit.e) * Time;

DeltaV_analytic = abs(final_orbit.omega - initial_orbit.omega) * ...
    sqrt(constants.mu_dim / initial_orbit.a) * 2 * initial_orbit.e / sqrt(1-initial_orbit.e^2) ...
    * 1/ (3 - 2 * initial_orbit.e^2) ;



% DeltaV_analytic2 = abs(final_orbit.omega - initial_orbit.omega) * 2 * initial_orbit.e * pi / ...
%     sqrt(1-initial_orbit.e^2) *  sqrt(constants.mu_dim / initial_orbit.a)  * ...
%     1/(-1 + initial_orbit.e + (60*pi-48*pi*initial_orbit.e^2)/(12*(1-initial_orbit.e^2))) 

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


% =========================================================================
% Analytic equations
% =========================================================================
figure
plot(Time/86400, omega_analytic*180/pi,'LineWidth',2)
xlabel('Time [days]')
ylabel('\omega [deg]')  
grid on
title('Analytic equation')
                     