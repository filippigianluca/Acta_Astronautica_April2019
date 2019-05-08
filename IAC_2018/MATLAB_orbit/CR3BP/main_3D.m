addpath(genpath('../spaceart_toolbox'))

clear; close all; clc;


%% Constants

% Earth gravitational parameter [km^3/s^2]
Earth.mu = astro_constants(13);

% Earth mass [kg]
Earth.mass = astro_constants(13)/astro_constants(1);

% Moon gravitational parameter [km^3/s^2]
Moon.mu = astro_constants(20);

% Moon mass [kg]
Moon.mass = astro_constants(20)/astro_constants(1);

% Moon distance from Earth [km]
Moon.mean_distance = 384388;

% Adimensional length and time for Earth-Moon rotating reference frame
L_star = Moon.mean_distance;
T_star = sqrt( L_star^3 / (astro_constants(1) * (Earth.mass + Moon.mass)));




%% CR3BP parameters
% Reference: Curtis, Orbital Mechanics for Engineering Students (2.12)

% Earth-Moon distance [km]
r12 = Moon.mean_distance;

% Earth and Moon masses [kg]
m1 = Earth.mass;
m2 = Moon.mass;

% Total mass of the system [kg]
M = m1 + m2;

% Dimensionless mass rations []
pi1 = m1 / M;
pi2 = m2 / M;

% Position of mass 1 in the rotated reference frame
x1 = -pi2 * r12;

% Position of mass 2 in the rotated reference frame
x2 =  pi1 * r12;


%% Numerical integration

% Adimensional initial conditions: x0, y0, z0, vx0, vy0, vz0 (IC for DRO)
% pos0_vel0 = [1.235730563909383, 0, 0, 0,  -0.556662961164044, 0];
pos0_vel0 = [0.1, 0, 0.01, 0, 3.91, 0];

% Adimensional time of flight
ToF =  2.071595961347605;

% Numerical integration
[x, y, z, vx, vy, vz, t] = CR3BPIntegration3D_adim(pos0_vel0(1,1), pos0_vel0(1,2), pos0_vel0(1,3), ...
                                                   pos0_vel0(1,4), pos0_vel0(1,5), pos0_vel0(1,6), ...
                                                   linspace(0, 5*ToF, 1001), Earth, Moon, L_star, T_star);


%% Plot
figure
plot3(x*L_star, y*L_star, z*L_star, 'k')
hold on
plot3(x1, 0, 0,'ko','MarkerFaceColor','k')
plot3(x2, 0, 0,'ko','MarkerFaceColor','k')
grid on
xlabel('x (ROTATING) [km]')
ylabel('y (ROTATING) [km]')
zlabel('z (ROTATING) [km]')
title('Rotating reference frame')
axis equal
cameratoolbar