addpath(genpath('../../spaceart_toolbox'))

% clear; close all; clc;

% reference orbit
rp = 6678; 
ra = 384400;
i = 7*pi/180;
Om = 0*pi/180;
om = 0*pi/180;
ni0 = 3*pi/4;

% initial (and constant except for M if da is not 0) difference for SC 1
da1 = 0.0;
de1 = 0.0;
di1 = -0.002*pi/180;
dOm1 = -0.0004*pi/180;
dom1 = 0.0*pi/180;
dM01 = 0.0*pi/180;

% initial (and constant except for M if da is not 0) difference for SC 2
da2 = 0.0;
de2 = 0.0;
di2 = 0.002*pi/180;
dOm2 = 0.0004*pi/180;
dom2 = 0.0*pi/180;
dM02 = 0.0*pi/180;

% computing other parameters
a = (ra+rp) / 2; 
e = (ra-rp) / (ra+rp);
mu = 398600.4418;
n = sqrt(mu / a^3);
T = 2 * pi / n;
dn1 = sqrt(mu / (a + da1)^3) - n;
dn2 = sqrt(mu / (a + da2)^3) - n;

E0 = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(ni0 / 2));
if(E0 < 0)
    E0 = E0 + 2 * pi;
end
M0 = E0 - e * sin(E0);
Dt0 = M0 / ni0;

thetas = linspace(ni0, pi, 1001); % vector of value for the reference true anomaly
nb = length(thetas);
rel_pos = zeros(3, nb);
posvel1 = rel_pos;
posvel2 = rel_pos;

time = zeros(size(thetas));

for k = 1 : nb
    nrev = floor(thetas(k) / (2 * pi));
    theta = thetas(k) - 2 * pi * nrev;
    E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(theta / 2));
    if(E < 0)
        E = E + 2 * pi;
    end
    M = E - e * sin(E);
    Dt = M / n + nrev * T;
    
    time(k) = Dt - Dt0;
    
    dM1 = dM01 + dn1 * Dt;
    posvel1(:, k) = lrom([a; e; i; Om; om; theta], [da1; de1; di1; dOm1; dom1; dM1]);
    
    dM2 = dM02 + dn2 * Dt;   
    posvel2(:, k) = lrom([a; e; i; Om; om; theta], [da2; de2; di2; dOm2; dom2; dM2]);
    
    rel_pos(:, k) = posvel1(1:3,k) - posvel2(1:3, k);
end

thetasbis = thetas * 180 / pi;
rel_dist = sqrt(rel_pos(1, :).^2 + rel_pos(2, :).^2 + rel_pos(3, :).^2);

figure; hold on; grid;
plot(thetasbis, rel_dist, 'b', 'LineWidth', 2);
axis([min(thetasbis) max(thetasbis) -Inf Inf]);
xlabel('true anomaly (deg)'); ylabel('relative distance (km)');

figure; hold on; grid;
plot(time, rel_dist, 'b', 'LineWidth', 2);
axis([min(time) max(time) -Inf Inf]);
xlabel('time (s)'); ylabel('relative distance (km)');

return

% plots of relative position in local orbital frame
figure; hold on;

subplot(3,1,1); hold on; grid;
plot(thetasbis, rel_pos(1, :), 'b', 'LineWidth', 2);
axis([min(thetasbis) max(thetasbis) -Inf Inf]);
xlabel('true anomaly (deg)'); ylabel('radial (km)');

subplot(3,1,2); hold on; grid;
plot(thetasbis, rel_pos(2, :), 'b', 'LineWidth', 2);
axis([min(thetasbis) max(thetasbis) -Inf Inf]);
xlabel('true anomaly (deg)'); ylabel('transversal (km)');

subplot(3,1,3); hold on; grid;
plot(thetasbis, rel_pos(3, :), 'b', 'LineWidth', 2);
axis([min(thetasbis) max(thetasbis) -Inf Inf]);
xlabel('true anomaly (deg)'); ylabel('out-of-plane (km)');

figure; hold on; grid;
plot(posvel1(2, :), posvel1(3, :), 'b', 'LineWidth', 2);
plot(posvel2(2, :), posvel2(3, :), 'r', 'LineWidth', 2);
axis equal;

figure; hold on; grid;
plot(rel_pos(2, :), rel_pos(3, :), 'k', 'LineWidth', 2);
