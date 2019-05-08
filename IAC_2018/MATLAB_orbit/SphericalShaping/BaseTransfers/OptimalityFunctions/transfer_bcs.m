function PSI = transfer_bcs(X0,Xf,r0,rf,v0,vf)
%
% 
%transfer_bcs: Boundary Condition Function for the Two Point Transfer.
% Function that writes the boundaries conditions for the two
% points transfer. They have to be written as residuals (expressions that
% go to zero)
%
% INPUT
% X0: vector with initial boundary conditions (r,v,lambdar,lambdav)
% Xf: vector with final boundary conditions (r,v,lambdar,lambdav)
% r0,rf,v0,vf: initial and final values of position and velocity to match

% OUTPUT
% PSI: vector with the boundary conditions residuals (variables - known
% values of boundary conditions). Column vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a column vector with the 8 boundary conditions on the state &
% costate variables

% Transform into column vector the boundary values
r0 = r0(:);
rf = rf(:);
v0 = v0(:);
vf = vf(:);

% Write PSI (residual)
PSI = [X0(1:3) - r0;
       X0(4:6) - v0;
       Xf(1:3) - rf;
       Xf(4:6) - vf];

return

