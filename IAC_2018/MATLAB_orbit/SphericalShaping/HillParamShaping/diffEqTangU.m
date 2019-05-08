function [xdot] = diffEqTangU(t,x,mu,uMax)

%diffEqTangU: function that writes the differential equations for
%the optimal control problem (minimum Delta V) using a tangential thrust
%equal in modulus to uMax
%
% INPUT
% t: time
% x: cartesian state (r,v) and Lagrange multipliers (r,v). Column vector
% mu: gravitational parameter
% uMax: modulus of tangential thrust

% OUTPUT
% xdot: time derivative of cartesian state (r,v) and of Lagrange multipliers (r,v). Column vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write the time derivatives

% Write the cartesian state and adjoints
r = x(1:3);
v = x(4:6);

% Find the unit vector correspondent to the velocity
modV = norm(v);
unitV = v/modV;

% Obtain the tangential control, of modulus uMax and direction V
uVec = uMax*unitV;

% Write the differential equations (from Homotopic approach and
% pseudospectral method...)
rdot = v;
vdot = -mu/norm(r)^3*r + uVec;

% Write the time derivatives as outputs
xdot(1:3) = rdot;
xdot(4:6) = vdot;

% Transpose for requirement of ode45 (and change to [1 1 1]'): output is a
% column vector
xdot = xdot(:);
end

