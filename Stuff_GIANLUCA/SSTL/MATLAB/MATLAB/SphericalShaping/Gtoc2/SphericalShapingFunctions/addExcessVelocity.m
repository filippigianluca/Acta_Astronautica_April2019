function [VinfCart] = addExcessVelocity(Vinf,alpha,beta,DepBody,Theta_dep)
% Function that adds an initial excess velocity to the departure from the
% planets, useful for some GTOC problems

% INPUT
% Vinf = modulus of excess velocity
% alpha = azimuth from tangential direction
% beta = elevation from orbit plane
% DepBody = structure with the Kepler elements of departure planet
% Theta_dep = true anomaly at departure

% OUTPUT
% VinfCart = vector with excess velocity in cartesian coordinates [X,Y,Z]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the velocity in tangential, normal, out of plane frame
Vt = Vinf*cos(beta)*cos(alpha);
Vn = Vinf*cos(beta)*sin(alpha);
Vh = Vinf*sin(beta);

% Calculate the flight path angle at departure
gamma = atan( (DepBody.e*sin(Theta_dep))/(1 + DepBody.e*cos(Theta_dep)) );

% Calculate the velocity in radial and transverse reference frame
Vr = Vt*sin(gamma) + Vn*cos(gamma);
Vtr = Vt*cos(gamma) - Vn*sin(gamma);

% Calculate the velocity in the perifocal frame
V_xi = Vr*cos(Theta_dep) - Vtr*sin(Theta_dep);
V_eta = Vr*sin(Theta_dep) + Vtr*cos(Theta_dep);

V_pf = [V_xi; V_eta; Vh];

% Transform the vector from perifocal frame to geocentric frame
[VinfCart] = KeplElem2rvMatrix(DepBody.Omega, DepBody.i, DepBody.w, V_pf);


end

