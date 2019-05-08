function sp = dc_eq(t,s,u,kp)

% Equation of the dynamics of the two body problem, with control.
%
%   sp = dc_eq(t,s,u,kp)
%
% INPUT:
%       t = time
%       s = state vector (position, velocity)
%       u = control vector
%       kp = G*mp con mp = mass of the attracting body
%                     G = universal gravity constant
%
% OUTPUT:
%       sp = derivative of the state
%
% FUNCTIONS CALLED:
%  (none)
%
% - Camilla Colombo - 07/03/2005
% Revised by Matteo Ceriotti - 10-01-2007
%
% ------------------------- - SpaceART Toolbox - --------------------------

r = s(1:3);
v = s(4:6);
u = u(:);

sp = [       v        ;...
      -kp/norm(r)^3*r + u];

return