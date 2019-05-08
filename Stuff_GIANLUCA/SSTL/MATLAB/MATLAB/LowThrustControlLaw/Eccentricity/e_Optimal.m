% =========================================================================
% Optimal thrusting angles strategy for change in the eccentricity
% Elevation equal to zero and azimuth in a specific direction
% Only numerical integration is realized in this file.
% References: 
% A. E. Petropolous, "Simple Control Laws for low-thrust orbit
% transfers", AAS/AIAA Astrodynamics Specialist Conference, 2003
% =========================================================================
% Marilena Di Carlo
% marilena.di-carlo@strath.ac.uk


clear
close all
addpath(genpath('..\..\spaceart_toolbox'))
addpath('../')


%% Change of e and i

mode = 'optimum_eccentricity';

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
% Input - user defined
% =========================================================================

% Spacecraft mass [kg]
m0 = 1000;

% Spacecraft thrust [N]
T = 0.08;

% Acceleration [m/s^2]
f = T / m0;

% Acceleration [km/s^2]
f = f * 1e-3;

% Specific impulse [s]
Isp = 3000;

% Specific impulse [TU]
% Isp = Isp / constants.TU;

% Initial orbit elements [km, rad]
initial_orbit.a     = 7578;
initial_orbit.e     = 0.00000001;
initial_orbit.i     = 30*pi/180;
initial_orbit.omega = 0*pi/180;
initial_orbit.Omega  = 0*pi/180;
initial_orbit.E = 10*pi/180;

% Final orbit [km, rad]
final_orbit.e = 0.1188;

if final_orbit.e < initial_orbit.e
    f = -f;
elseif final_orbit.e == initial_orbit.e
     warning('Initial and final eccentricity are equal.');
end


element_target = 'eccentricity';
beta = 0;
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.e));
tic
[Time, Kep_Elements] = ode45(@(t,x)rate_of_change_E(t,x,f,beta,constants, mode),[0 365*86400], ...
                        [initial_orbit.a initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega initial_orbit.E], ...
                        options_ODE);
                    
toc
     
Delta_V_numeric = Time(end) * f;

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
                    

