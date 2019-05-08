% Test the verification algorithm: generate a random shape (a spiral) and
% compute with finite differences its velocity and acceleration. Then
% compute the corresponding control acceleration. Then integrate the
% dynamics with that control acceleration.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation
clc
clear all
close all

d2r = pi/180;
r2d = 180/pi;
AU2km = 149597870.66;                         % [km]
km2AU = 1/AU2km;

% Set constants
T_sid = 86164.1004;                           % [s]
mu_S = 1.327178*1e11*km2AU^3*T_sid^2;         % [AU^3/day^2]

% Spiral
a = 1;
b = 0.5;
theta = [0:pi/200:5*pi];

time = linspace(1,1000,length(theta));

r = a+b*theta;

% Position
Rx = r.*cos(theta);
Ry = r.*sin(theta);
Rz = 0*Rx;

% Velocity
Vx = diff(Rx)./diff(time);
Vx(end+1) = 2*Vx(end)-Vx(end-1);

Vy = diff(Ry)./diff(time);
Vy(end+1) = 2*Vy(end)-Vy(end-1);

Vz = 0*Vx;

% Acceleration
Ax = diff(Vx)./diff(time);
Ax(end+1) = 2*Ax(end)-Ax(end-1);

Ay = diff(Vy)./diff(time);
Ay(end+1) = 2*Ay(end)-Ay(end-1);

Az = 0*Ax;

% Plot
figure;
plot3(Rx, Ry, Rz,'g','LineWidth',2.0);

figure;
subplot(3,1,1);
plot(time,Rx,'b');
subplot(3,1,2);
plot(time,Vx,'g');
subplot(3,1,3);
plot(time,Ax,'r');
suptitle('x');

figure;
subplot(3,1,1);
plot(time,Ry,'b');
subplot(3,1,2);
plot(time,Vy,'g');
subplot(3,1,3);
plot(time,Ay,'r');
suptitle('y');

R = [Rx' Ry' Rz'];
RMod = (Rx.^2 + Ry.^2 + Rz.^2)';
V = [Vx' Vy' Vz'];
A = [Ax' Ay' Az'];

% Control
for i=1:3
   U(:,i) = A(:,i) + (mu_S./RMod.^3).*R(:,i);
end

figure;
subplot(3,1,1);
plot(time,A(:,1),'b');
subplot(3,1,2);
plot(time,U(:,1),'g');
subplot(3,1,3);
plot(time,(mu_S./RMod.^3).*R(:,1),'r');
suptitle('Acc x');

figure;
subplot(3,1,1);
plot(time,A(:,2),'b');
subplot(3,1,2);
plot(time,U(:,2),'g');
subplot(3,1,3);
plot(time,(mu_S./RMod.^3).*R(:,2),'r');
suptitle('Acc y');

% Verification
r0 = R(1,:);
v0 = V(1,:);
[R_ver,V_ver] = TrajectoryFromControl(time,r0,v0,U',mu_S,R',V');