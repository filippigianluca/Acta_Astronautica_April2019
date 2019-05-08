addpath(genpath('../spaceart_toolbox'))

clear; close all; clc;


%% Constants

% Earth gravitational parameter [km^3/s^2]
Earth.mu = astro_constants(13);

% Earth mass [kg]
Earth.mass = astro_constants(13)/astro_constants(1);

% Earth radius [km]
Earth.radius = astro_constants(23);

% Moon gravitational parameter [km^3/s^2]
Moon.mu = astro_constants(20);

% Moon mass [kg]
Moon.mass = astro_constants(20)/astro_constants(1);

% Moon radius [km]
Moon.radius = astro_constants(30);

% Moon distance from Earth [km]
Moon.mean_distance = 384388;

% Adimensional length and time for Earth-Moon rotating reference frame
L_star = Moon.mean_distance;
T_star = sqrt( L_star^3 / (astro_constants(1) * (Earth.mass + Moon.mass)));

% Coordinates of a sphere
[Xs, Ys, Zs] = sphere(20);

% True anomalies
theta = linspace(0*pi/180, 360*pi/180, 1001);



%% CR3BP parameters
% Reference: Curtis, Orbital Mechanics for Engineering Students (2.12)

% Earth-Moon distance [km]
r12 = Moon.mean_distance;

% Earth and Moon masses [kg]
m1 = Earth.mass;
m2 = Moon.mass;

% Total mass of the system [kg]
M = m1 + m2;

% Dimensionless mass ratios []
pi1 = m1 / M;
pi2 = m2 / M;

% Position of mass 1 in the rotated reference frame
x1 = -pi2 * r12;

% Position of mass 2 in the rotated reference frame
x2 =  pi1 * r12;


%% Initial conditions

% Orbital parameters of the Moon in ECI
kep_moon = [4e5 0.0 27*pi/180 0*pi/180 0*pi/180 0*pi/180];

% Orbital period of the system
T_sys = 2*pi*sqrt(kep_moon(1)^3/Earth.mu);

% Rotation from ECI to rotating frame
R_ni0  = [cos(kep_moon(6)) -sin(kep_moon(6)) 0; sin(kep_moon(6)) cos(kep_moon(6)) 0; 0 0 1];
R_argp = [cos(kep_moon(5)) -sin(kep_moon(5)) 0; sin(kep_moon(5)) cos(kep_moon(5)) 0; 0 0 1];
R_inc  = [1 0 0; 0 cos(kep_moon(3)) -sin(kep_moon(3)); 0 sin(kep_moon(3)) cos(kep_moon(3))];
R_raan = [cos(kep_moon(4)) -sin(kep_moon(4)) 0; sin(kep_moon(4)) cos(kep_moon(4)) 0; 0 0 1];
ECI2ROT = (R_raan*R_inc*R_argp*R_ni0).';

% Initial conditions of the spacecraft in the ECI frame
ra_parking   = Earth.radius+500; %((86400/(2*pi))^2*Earth.mu)^(1/3); %
rp_parking   = Earth.radius+500;
sma_parking  = (ra_parking+rp_parking)/2;
ecc_parking  = 1-rp_parking/sma_parking;
inc_parking  = 0*pi/180;
raan_parking = 0*pi/180;
argp_parking = 0*pi/180;
ni0_parking  = 0*pi/180;
T_parking = 2*pi*sqrt(sma_parking^3/Earth.mu);
kep_parking  = [sma_parking ecc_parking inc_parking raan_parking argp_parking ni0_parking];

X0_ECI = kep2cart(kep_parking, Earth.mu);

% Rotation to rotating frame
X0_rot(1:3,1) = ECI2ROT*X0_ECI(1:3).';
X0_rot(4:6,1) = ECI2ROT*X0_ECI(4:6).';

% Translation to the position of the Earth
X0_rot(1) = X0_rot(1)+x1;

% Removing centrifugal velocity (+ or - ?)
X0_rot(4:6) = X0_rot(4:6) -cross([0 0 2*pi/T_sys], X0_rot(1:3)).';

% Non-dimensionalisation
X0_rot = X0_rot/L_star;
X0_rot(4:6) = X0_rot(4:6)*T_star;

% Adimensional initial conditions: x0, y0, z0, vx0, vy0, vz0 (IC for DRO)
% pos0_vel0 = [1.235730563909383, 0, 0, 0,  -0.556662961164044, 0];
% pos0_vel0 = [0.2, 0, 0.01, 0, 2.5, 0];

pos0_vel0 = X0_rot;



%% Control sequence

control_sequence(1).time = 0.01*T_parking/T_star;
control_sequence(1).deltaV = [0, 3.2, 0]*T_star/L_star;

control_sequence(2).time = 2*T_parking/T_star;
control_sequence(2).deltaV = [0, 0, 0]*T_star/L_star;

control_sequence(3).time = 3*T_parking/T_star;
control_sequence(3).deltaV = [0, 0, 0]*T_star/L_star;

ctrl = reshape([vertcat(control_sequence.time) vertcat(control_sequence.deltaV)].', 1, []);


%% Numerical integration

% Adimensional time of flight
ToF =  2.071595961347605;

% Numerical integration
[x, y, z, vx, vy, vz, t] = CR3BPIntegration3D_adim_imp_ctrl(pos0_vel0(1), pos0_vel0(2), pos0_vel0(3), ...
                                                            pos0_vel0(4), pos0_vel0(5), pos0_vel0(6), ...
                                                            ctrl,...
                                                            [0, 6*ToF], Earth, Moon, L_star, T_star);


%% Plot
figure
plot3(x*L_star, y*L_star, z*L_star, 'k')
hold on
% plot3(x1, 0, 0,'ko','MarkerFaceColor','k')
surf(Xs*Earth.radius+x1, Ys*Earth.radius, Zs*Earth.radius, 'EdgeAlpha', 0)
% plot3(x2, 0, 0,'ko','MarkerFaceColor','k')
surf(Xs*Moon.radius+x2, Ys*Moon.radius, Zs*Moon.radius, 'EdgeAlpha', 0)
grid on
xlabel('x (ROTATING) [km]')
ylabel('y (ROTATING) [km]')
zlabel('z (ROTATING) [km]')
title('Rotating reference frame')
axis equal
cameratoolbar
