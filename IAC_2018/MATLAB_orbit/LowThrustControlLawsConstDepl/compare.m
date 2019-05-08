%==========================================================================

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
ToF_total = (550:10:1000)*86400;

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

for i = 1 : length(ToF_total)

% Change this: now two files might be called, the second called when in the
% first no solution to the problem is found. This should be avoided and
% everything should be in the same file
[DeltaV1(i), x_alpha, beta, ToF_1st, ToF_2nd] = RAANJ2_ai_2arcs(initial_orbit, final_orbit, ...
                                    ToF_total(i), ...
                                    f, constants, plot_flag);
                                
                                                    %                                 
if isnan(DeltaV1(i))
    warning('Transfer has no solution for given initial and final orbital elements and given time of flight')
end

f_direction = 1;
[best_solution, all_solutions] = ai_RAAN_2arcs(initial_orbit, final_orbit, ...
                                    ToF_total(i), f_direction, ...
                                    f, constants);
       
% -------------------------------------------------------------------------
% Minimum DeltaV solution
% -------------------------------------------------------------------------
% DeltaV minimum DeltaV solution
DeltaV2(i)    = best_solution.DeltaV;


[DeltaV3(i), alpha1, ToF1, alpha2, ToF2] = aRAAN_i_2arcs(initial_orbit, final_orbit,...
    ToF_total(i), ...
                                    f, constants);


end

ls = 16;

figure
plot(ToF_total/86400,DeltaV1,'LineWidth',2)
hold on
plot(ToF_total/86400,DeltaV2,'LineWidth',2)
plot(ToF_total/86400,DeltaV3,'LineWidth',2)
grid on
xlabel('ToF_{total} [days]','FontSize',ls)
legend('First strategy','Second strategy','Third strategy')
ylabel('\Delta V [km/s]','FontSize',ls)
set(gca,'FontSize',ls)