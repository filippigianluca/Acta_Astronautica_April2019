function [xdot] = DiffHillEquationsWithControl(u,x,mu,F_r,F_u,F_n)

%diffEqWithControl: function that writes the differential equations for a
% lt spacecraft with a given control
%
% INPUT
% t: time
% x: [r; v_r; G; u; i; h];
% mu: gravitational parameter
% u: matrix with the control acceleration

% OUTPUT
% xdot: time derivative of cartesian state (r,v). Column vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Column vector

xdot = zeros(6,1);

% Write the differential equations
xdot(1) = x(2)/(x(3)/x(1)^2  - (x(1)*cot(x(5))*sin(u)/x(3)) * F_n); 
xdot(2) = (x(3)^2/x(1)^3 - mu/x(1)^2 + F_r) / (x(3)/x(1)^2 - (x(1)*cot(x(5))*sin(u)/x(3)) * F_n) ;
xdot(3) = (x(1) * F_u) / (x(3)/x(1)^2 - (x(1)*cot(x(5))*sin(u)/x(3)) * F_n);
xdot(4) = (x(3)/x(1)^2 - (x(1)*cot(x(5))*sin(u)/x(3)) * F_n) / (x(3)/x(1)^2 - (x(1)*cot(x(5))*sin(u)/x(3)) * F_n);
xdot(5) = ((x(1)*cos(u)/x(3)) * F_n) / (x(3)/x(1)^2 - (x(1)*cot(x(5))*sin(u)/x(3)) * F_n);
xdot(6) = ((x(1)*sin(u)/x(3)/sin(x(5))) * F_n) / (x(3)/x(1)^2 - (x(1)*cot(x(5))*sin(u)/x(3)) * F_n);
end