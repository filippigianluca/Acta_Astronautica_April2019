% =========================================================================
% Pollard control law for variation of omega without variation of other
% orbital elements (semimajor axis and eccentricity).
% DeltaV computed analytically. 
% Orbital elements integrated numerically.
% Analytical integration could be realised in FABLE. Need to compute new
% analytic integrals.
% Reference: Pollard, "Simplified Analysis of low0thrust orbital manuevers"
% ==========================================================================
% Marilena Di Carlo, 2016
% marilena.di-carlo@strath.ac.uk


clear
close all
addpath(genpath('..\..\spaceart_toolbox'))



%% Change omega

mode = 'omega';

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
% Input, defined by the user
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
initial_orbit.E = 180*pi/180;

% Final orbit [km, rad]
final_orbit.omega = 40 * pi/180;


% =========================================================================
% DeltaV and ToF computation
% =========================================================================

% Transfer for change in eccentricity and inclination
[Delta_V_equation, ToF] = Pollard_omega(initial_orbit, final_orbit, f, constants);



% =========================================================================
% Numerical integration for variation of orbital elements
% =========================================================================
if final_orbit.omega < initial_orbit.omega
    f = -f;
end

beta=0;
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12);
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 ToF], ...
                        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);
                    
Delta_V_numeric = Time(end) * abs(f);


% =========================================================================
% Analytic expression variation omeg
% =========================================================================
omega_analytic = initial_orbit.omega +  ...
    3 * sqrt(initial_orbit.a * (1 - initial_orbit.e^2) / constants.mu_dim) * ...
    f * Time / (2 * initial_orbit.e);
                    
% =========================================================================
% Plot
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
                     