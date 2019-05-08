%==========================================================================
% Variation of semimajor axis, inclination and right ascension with J2 in a
% given time of flight using low-thrust - First strategy
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
i_final =56 * pi/180;

% Initial RAAN
RAAN_initial = 0*pi/180;
% Final RAAN
RAAN_final   =150*pi/180;

% Total ToF [s]
ToF_total = 600*86400;

% Low thrust acceleration [km/s^2]
f = 1.5e-7;

% Plot?
plot_flag = 0;



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

% Change this: now two files might be called, the second called when in the
% first no solution to the problem is found. This should be avoided and
% everything should be in the same file
[DeltaV, x_alpha, beta, ToF_1st, ToF_2nd] = RAANJ2_ai_2arcs(initial_orbit, final_orbit, ...
                                    ToF_total, ...
                                    f, constants, plot_flag);
                                
                                                    %                                 
if isnan(DeltaV)
    warning('Transfer has no solution for given initial and final orbital elements and given time of flight')
end



%% Plot

% Vector time for the first phase
t = linspace(0, ToF_1st, 1000);

% During the first phase the RAAN changes only because of J2
for i = 1 : length(t)
    RAAN_1st(i) = RAAN_initial -1.5 * sqrt(constants.mu) * constants.J2 * constants.R_Earth^2 * ...
        cos(i_initial) * a_initial^(-7/2) * t(i);
end

% Vector time for the second phase
t2 = linspace(0, ToF_2nd, 1000);

% Velocity variation during the second phase
V_t = sqrt(constants.mu/a_initial) - 2 * f * cos(beta) * x_alpha * t2 / pi;

% Semimajor axis variation during the first phase (circular orbits)
a_t = constants.mu./V_t.^2;

% Inclination variation during the second phase
i_t = i_initial + tan(beta) * sin(x_alpha) / x_alpha * log(a_t./a_initial)./2;

% Parameters for the computation of the varition of RAAN during the second
% phase
k2 = - 4 * log(a_final/a_initial) / (i_final - i_initial);

k3 = exp(k2 * i_t) .* (k2 * cos(i_t) + sin(i_t)) - ...
    exp(k2 * i_initial) * (k2 * cos(i_initial) + sin(i_initial));

k1 = -3/4 * constants.mu * pi * constants.J2 * constants.R_Earth^2 / ...
    (a_initial^4 * f * sin(beta) * sin(x_alpha)) * ...
    exp(4 * log(a_final/a_initial) * i_initial / (i_final - i_initial));

% Variation of RAAN during the second phase
RAAN_2nd = mod(RAAN_1st(end) + k1 / (1+k2^2) * k3, 2*pi);


% -------------------------------------------------------------------------
% Semimajor axis second phase
% -------------------------------------------------------------------------
ls = 16;
figure(1)
hold on
plot(t/86400, a_initial*ones(1,length(t)), 'LineWidth',2)
plot(t(end)/86400 + t2/86400, a_t, 'LineWidth',2)
grid on
xlabel('Time [days]','FontSize',ls)
ylabel('Semimajor axis [km]','FontSize',ls)
set(gca,'fontsize',ls)

% -------------------------------------------------------------------------
% Inclination second phase
% -------------------------------------------------------------------------
figure(2)
hold on
plot(t/86400, i_initial*ones(1,length(t))*180/pi, 'LineWidth',2)
plot(t(end)/86400 + t2/86400, i_t*180/pi, 'LineWidth',2)
grid on
ylim([50 56])
xlabel('Time [days]','FontSize',ls)
ylabel('Inclination [deg]','FontSize',ls)
set(gca,'fontsize',ls)

% -------------------------------------------------------------------------
% RAAN second phase
% ------------------------------------------------------------------------
figure(3)
hold on
plot(t/86400, mod(RAAN_1st, 2*pi)*180/pi, 'LineWidth',2)
plot(t(end)/86400 + t2/86400, mod(RAAN_2nd, 2*pi)*180/pi, 'LineWidth',2)
grid on
xlabel('Time [days]','FontSize',ls)
ylabel('RAAN [deg]','FontSize',ls)
set(gca,'fontsize',ls)
    


