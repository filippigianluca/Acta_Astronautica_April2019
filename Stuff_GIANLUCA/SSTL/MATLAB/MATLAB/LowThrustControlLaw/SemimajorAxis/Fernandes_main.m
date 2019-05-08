% =========================================================================
% Thrusting law for the low-thrust variation of semimajor axis,
% eccentricity and argument of the perigee.
% User has to define the input (initial and final orbital elements,
% low-thrst acceleration).
% The script plot the variation of the orbital elements and return the
% DeltaV required (workspace)
% =========================================================================
% Reference: 
% da Silva Fernandes

% Marilena Di Carlo, 2015
% marilena.di-carlo@strath.ac.uk

clear
close all


% =========================================================================
% Constants
% =========================================================================

% Earth gravitational parameter [km^3/s^2]
constants.mu = 398600;

% Gravity acceleration [km/s^2]
constants.g0 = 9.8 * 1e-3;

%% ========================================================================
% INPUT
% =========================================================================
%% Engine and spacecraft

% Spacecraft mass [kg]
m0 = 2000;

% Power [W = kg m^2 / s^3]
engine.P_max = 8000;

% Power [kg km^2 /s^3]
engine.P_max = engine.P_max * 1e-6;

% Specific impulse [s]
engine.Isp = 3000;


%% Initial conditions

% Semimajor axis [km]
initial_orbit.a = 7000;

% Eccentricity
initial_orbit.e = 0.1;

% Perigee argument [rad]
initial_orbit.omega = 0 * pi/180;


%% Final conditions

% Semimajor axis [km]
final_orbit.a = 10000;

% Eccentricity
final_orbit.e = 1;

% Perigee argument [rad]
final_orbit.omega = 0 * pi/180;


%% Transfer time

% Time [days]
T = 45;

% Number of element in time vector
n = 100;


%% ========================================================================
% Compute transfer
% =========================================================================

% Transfer time [s]
T = T * 86400;
t = linspace(0, T, 100);

% DeltaV, variation of orbital elements, mass and acceleration during the
% transfer
[DeltaV, a, e, omega, m, acc, thrust] = Fernandes_a_e_omega(initial_orbit, final_orbit, t, engine, m0, constants);
 
% DeltaV =integral(@(tt)deltaV_Fernandes(tt, a0, k0, h0, pa0, pk0, ph0, E, C, mu,P_max,Isp,g0,m0), 0, t(end));
% DeltaV = Isp * g0 * log(m0 / m(end));
 
 

% =========================================================================
% Plot results
% =========================================================================

figure
subplot(2,2,1)
hold on
plot(t/86400,a,'LineWidth',2)
line([0 t(end)/86400],[final_orbit.a final_orbit.a],'Color','r','LineWidth',2)
grid on
xlabel('Time [days]')
ylabel('Semimajor axis [km]')
subplot(2,2,2)
hold on
line([0 t(end)/86400],[final_orbit.e final_orbit.e],'Color','r','LineWidth',2)
plot(t/86400,e,'LineWidth',2)
grid on
xlabel('Time [days]')
ylabel('Eccentricity')
subplot(2,2,3)
hold on
line([0 t(end)/86400],[final_orbit.omega*180/pi final_orbit.omega*180/pi],'Color','r','LineWidth',2)
plot(t/86400,omega*180/pi,'LineWidth',2)
grid on
xlabel('Time [days]')
ylabel('Perigee argument [deg]')
subplot(2,2,4)
plot(t/86400, m,'LineWidth',2)
grid on
xlabel('Time [days]')
ylabel('Mass [kg]')



figure
plot(t/86400, acc*1000, 'LineWidth',2)
xlabel('Time [days]')
ylabel('Acceleration [m/s^2]')


 figure
 plot(t/86400,thrust, 'LineWidth',2)
 xlabel('Time [days]')
 ylabel('Thrust [N]')

