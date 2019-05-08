% =========================================================================
% Pollard control law for simultaneous variation of eccentricity and
% inclination.
% DeltaV computed analytically. 
% Orbital elements integrated numerically.
% Reference: Pollard, "Simplified Analysis of low0thrust orbital manuevers"
% ==========================================================================
% Marilena Di Carlo
% marilena.di-carlo@strath.ac.uk


clear
close all
addpath(genpath('..\..\spaceart_toolbox'))
addpath(genpath('..\..\Propagator'))
addpath('..\')


%% Change of e and i

% Mode of thrusting
mode = 'eccentricity_inclination';

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
m0 = 2000;

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
initial_orbit.e     = 0.005;
initial_orbit.i     = 30*pi/180;
initial_orbit.omega = 0*pi/180;
initial_orbit.Omega  = 0*pi/180;
initial_orbit.E = 180*pi/180;

% Final orbit [km, rad]
final_orbit.e = 0.05;
final_orbit.i = 20*pi/180;



% =========================================================================
% DeltaV and ToF computation
% =========================================================================
% Transfer for change in eccentricity and inclination
[Delta_V_equation, ToF, beta] = Pollard_e_i(initial_orbit, final_orbit, f, constants);


% keyboard

% =========================================================================
% Numerical integration for variation of orbital elements
% =========================================================================

if final_orbit.e < initial_orbit.e && final_orbit.i < initial_orbit.i
    f = -f; 
elseif final_orbit.e > initial_orbit.e && final_orbit.i < initial_orbit.i
%     beta = -beta;
elseif final_orbit.e < initial_orbit.e && final_orbit.i == initial_orbit.i
    f = -f;
elseif final_orbit.e < initial_orbit.e && final_orbit.i > initial_orbit.i
    f = -f;
%     beta = -beta;
elseif final_orbit.e == initial_orbit.e
    warning('This transfer can not work with initial eccentricity equal to final eccentricity')
    return
end

options_ODE = odeset('AbsTol',1e-8,'RelTol',1e-8);

[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 ToF], ...
                        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);

Delta_V_numeric = Time(end) * abs(f);  


% =========================================================================
% Analytic equation for the variation of e
% =========================================================================
e_analytic =  sin( asin(initial_orbit.e) + 1.5 *cos(beta) * sqrt(initial_orbit.a / constants.mu_dim) * ...
    f * Time);
tic


% Analytic equations
num = cos(0.25 * (2 * asin(initial_orbit.e) + 3 * cos(beta) * f * ...
    sqrt(initial_orbit.a / constants.mu_dim) * Time )) - ...
    sin(0.25 * (2 * asin(initial_orbit.e) + 3 * cos(beta) * f * ...
    sqrt(initial_orbit.a / constants.mu_dim) * Time ));

den = cos(0.25 * (2 * asin(initial_orbit.e) + 3 * cos(beta) * f * ...
    sqrt(initial_orbit.a / constants.mu_dim) * Time )) + ...
    sin(0.25 * (2 * asin(initial_orbit.e) + 3 * cos(beta) * f * ...
    sqrt(initial_orbit.a / constants.mu_dim) * Time ));

num0 = cos( asin(initial_orbit.e/2) ) - sin( asin(initial_orbit.e/2) );

den0 =  cos( asin(initial_orbit.e/2) ) + sin( asin(initial_orbit.e/2) );

i_analytic = initial_orbit.i - ...
    4/3 * cos(initial_orbit.omega) * sin(beta) /(pi * cos(beta)) * ...
    ( 2 * log(num./den)  + ...
    sin(initial_orbit.e + 1.5 * cos(beta) * f * Time * sqrt(initial_orbit.a / constants.mu_dim))...
    ) + ...
    + 4/3 * cos(initial_orbit.omega) * sin(beta) / (pi * cos(beta)) * ...
    (2 * log(num0 / den0) + initial_orbit.e );


% =========================================================================
% Plot numeric integration
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
% Plot analytic equations
% =========================================================================
figure
subplot(1,2,1)
plot(Time/86400, e_analytic,'LineWidth',2)
grid on
xlabel('Time [days]')
ylabel('e')  
title('Analytic equation')
subplot(1,2,2)
plot(Time/86400, i_analytic*180/pi,'LineWidth',2)
grid on
xlabel('Time [days]')
ylabel('i [deg]')  
title('Analytic equation')
