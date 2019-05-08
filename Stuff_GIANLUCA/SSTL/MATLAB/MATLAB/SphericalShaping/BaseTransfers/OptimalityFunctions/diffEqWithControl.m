function [xdot] = diffEqWithControl(t,x,mu,u,tVec)

%diffEqWithControl: function that writes the differential equations for
%the optimal control problem (minimum Delta V) with the control given as
%input.
%
% INPUT
% t: time
% x: cartesian state (r,v). Column vector
% u: matrix with control
% mu: gravitational parameter
% OUTPUT
% xdot: time derivative of cartesian state (r,v) and of Lagrange multipliers (r,v). Column vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write the time derivatives

% Interpolate to obtain u at time t of integration, from known control
% vector u that comes from the direct method solution
U(1) = spline(tVec,u(1,:),t);
U(2) = spline(tVec,u(2,:),t);
U(3) = spline(tVec,u(3,:),t);
  
% Write the cartesian state
r = x(1:3);
v = x(4:6);

% Write the differential equations
rdot = v;
vdot = -mu/norm(r)^3*r + U';

% Write the time derivatives as outputs
xdot(1:3) = rdot;
xdot(4:6) = vdot;

% Transpose for requirement of ode45, to have a column vector
xdot = xdot(:);
end

