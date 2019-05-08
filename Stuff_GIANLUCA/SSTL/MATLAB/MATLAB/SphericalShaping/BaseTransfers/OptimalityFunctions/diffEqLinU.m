function [xdot] = diffEqLinU(t,x,mu,uMax)

%diffEq: function that writes the differential equations for
%the optimal control problem (minimum Delta V)
%
% INPUT
% t: time
% x: cartesian state (r,v) and Lagrange multipliers (r,v). Column vector
% mu: gravitational parameter

% OUTPUT
% xdot: time derivative of cartesian state (r,v) and of Lagrange multipliers (r,v). Column vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write the time derivatives

% Write the cartesian state and adjoints
r = x(1:3);
v = x(4:6);
lambdar = x(7:9);
lambdav = x(10:12);

% Find the direction of the control from the adjoints
uDir = - lambdav/norm(lambdav);

% Calculate the Switching Function
% Param Isp g0 m
T_sid = 86164.1004; 
AU2km = 149597870.66;  
Ispgom = (3000*9.81*1e-3)*AU2km/T_sid;
S = 1 - Ispgom*norm(lambdav);

% Choose control modulus depending upon the Switching Function
if S < 0 
    uMod = uMax;
else
    uMod = 0;
end

% Define the vector of control
uVec = uMod*uDir;

% % Write the differential equations (from Vasile paper)
% rdot = v;
% vdot = -mu/norm(r)^3*r + uVec;
% lambdardot = -lambdav.*(3*mu/norm(r)^5*(r.^2) - mu/norm(r)^3*[1 1 1]');
% lambdavdot = -lambdar;

% Write the differential equations (from Homotopic approach and
% pseudospectral method...)
rdot = v;
vdot = -mu/norm(r)^3*r + uVec;
lambdardot = lambdav*mu/norm(r)^3 - 3*mu/norm(r)^5*(dot(lambdav,r))*r;
lambdavdot = -lambdar;

% Write the time derivatives as outputs
xdot(1:3) = rdot;
xdot(4:6) = vdot;
xdot(7:9) = lambdardot;
xdot(10:12) = lambdavdot;

% Transpose for requirement of ode45 (and change to [1 1 1]'): output is a
% column vector
xdot = xdot(:);
end

