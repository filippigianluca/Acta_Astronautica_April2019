% =========================================================================
% Optimal thrusting angles strategy for change in the argument of perigee.
% Reference: Simple Control Laws for LowThrust Orbit Transfers, Petropolous
% Numerical integration of the orbital elements.
% Analytical integration could be realised with FABLE. New analytical
% integrals to be computed.
% =========================================================================
% Marilena Di Carlo, 2015
% marilena.di-carlo@strath.ac.uk



clear
close all
addpath(genpath('spaceart_toolbox'))
addpath('../')


% Mode for the thrust control
mode = 'optimum_omega';


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
initial_orbit.omega = 00*pi/180;
initial_orbit.Omega  = 30*pi/180;
initial_orbit.E = 10*pi/180;

% Final orbit [km, rad]
final_orbit.omega = 20*pi/180;

% Time for the integration - the integration stops when the desired element
% value is reached
time_int = 365 * 86400;

% =========================================================================
% Numerical Integration
% =========================================================================

% An event function is used to stop the integration when the targeted
% element has reached the desired value

if final_orbit.omega < initial_orbit.omega
    f = -f;
end
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.omega));
tic
beta = 0;
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 time_int], ...
                        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);
                  
Delta_V_numeric = Time(end) * abs(f);


% Delta_V_analytic = abs(final_orbit.omega - initial_orbit.omega) / ...
%     (0.5 * sqrt(initial_orbit.a / constants.mu) * sqrt(1-initial_orbit.e^2) * 1/initial_orbit.e^4 * ...
%     (4 - 2 * initial_orbit.e^2 + 9 *initial_orbit.e^4));
% 
% k = -1+initial_orbit.e + (60 * pi - 48 * pi * initial_orbit.e^2)/(12*(1-initial_orbit.e^2))
% DeltaV_analytic = abs(final_orbit.omega - initial_orbit.omega) *1000 / ...
%     (1/(2*pi) * sqrt(1-initial_orbit.e^2)/initial_orbit.e * sqrt(initial_orbit.a/constants.mu) * k)
% 
% DeltaV_analytic = abs(final_orbit.omega - initial_orbit.omega) * ...
%     sqrt(initial_orbit.a/constants.mu_dim) * ...
%     2 * initial_orbit.e / (sqrt(1-initial_orbit.e^2) * (3-2*initial_orbit.e^2));
% 
% fe = sqrt(1-initial_orbit.e^2)/initial_orbit.e * ...
%     (4*pi + 2*pi/(initial_orbit.e^2) + 3/2 + initial_orbit.e/3 + ...
%     (1+initial_orbit.e^2)/initial_orbit.e - ...
%     (1-initial_orbit.e^2)/initial_orbit.e^2 * log(1+initial_orbit.e));
% 
% DeltaV = abs(final_orbit.omega - initial_orbit.omega) * 2 * pi / ...
%     initial_orbit.a^2 * fe;
% 
% par = 1.5 + initial_orbit.e/3 + (1+initial_orbit.e^2)/initial_orbit.e + ...
%     2*pi/initial_orbit.e^2 + 4 * pi - (1-initial_orbit.e^2) / initial_orbit.e^2 * ...
%     log(1 + initial_orbit.e);
% 
% p = initial_orbit.a * (1-initial_orbit.e^2);
% h = sqrt(p * constants.mu_dim);
% T = 2 * pi * sqrt(initial_orbit.a^3/ constants.mu_dim);
%  abs(final_orbit.omega - initial_orbit.omega) * h^2 * initial_orbit.e^3 * T / ...
%      (4*pi*p^3)

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
                     