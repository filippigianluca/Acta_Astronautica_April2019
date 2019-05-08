function [xdot] = diffEqScaled(t,x,mu,rS,vS,lrS,lvS,tfS)

%diffEq: function that writes the differential equations for
%the optimal control problem (minimum Delta V). Version for the scaled
%equations
%
% INPUT
% t: time
% x: cartesian state (r,v) and Lagrange multipliers (r,v). Column vector
% mu: gravitational parameter
% rS,vS,lrS,lvS,tfS: scaling parameters (the scaled variables are r_tilde =
% r/rS)
%
%
%
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

% % Write the differential equations (from Vasile paper)
% rdot = v;
% vdot = -mu/norm(r)^3*r - lambdav;
% lambdardot = -lambdav.*(3*mu/norm(r)^5*(r.^2) - mu/norm(r)^3*[1 1 1]');
% lambdavdot = -lambdar;

% Write the differential equations (from Homotopic approach and
% pseudospectral method...)
rdot = v*(tfS*vS/rS);
vdot = -mu/norm(r)^3*r*(tfS/vS/rS^2) - lambdav*(tfS*lvS/vS);
lambdardot = (lambdav*mu/norm(r)^3 - 3*mu/norm(r)^5*(dot(lambdav,r))*r)*(tfS*lvS/lrS/rS^3);
lambdavdot = -lambdar*(tfS*lrS/lvS);

% Write the time derivatives as outputs
xdot(1:3) = rdot;
xdot(4:6) = vdot;
xdot(7:9) = lambdardot;
xdot(10:12) = lambdavdot;

% Transpose for requirement of ode45 (and change to [1 1 1]'): output is a
% column vector
xdot = xdot(:);
end

