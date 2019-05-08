function rho = lrom(kep_rob,doe)

% Linearized Relative Orbit Motion for general elliptic orbits
%
% rho = lrom(kep_rob,doe)
%
% From Schaub p. 619.
% It calculates the relative position coordinates in terms of the orbit
% element differences. They are expressed in a radial-transverse-h
% {r,theta,h} reference frame.
% All the angles in [rad], other dimensions to be consistent.
%
% INPUT:
%       kep_rob = orbital parameters of the chief orbit [a e i Om om f]
%                 f is the true anomaly where you want to know the relative
%                 position
%       doe = orbit element differences [da de di dOm dom dM]
%             dM is the difference in mean anomaly.
%
% OUTPUT:
%       rho = relative position vector expressed in {r,theta,h} reference
%             frame:
%               r-axis: direction of the orbit radius
%               h-axis: direction of angular momentum
%               theta-axis: transverse direction in the orbit plane,
%                           completes the reference frame
%
% functions called: none
%
% - Camilla Colombo - 07/06/2006
% - Revised by Matteo Ceriotti - 19/05/2007
%
% ------------------------- - SpaceART Toolbox - --------------------------

a  = kep_rob(1);     da  = doe(1);
e  = kep_rob(2);     de  = doe(2);
i  = kep_rob(3);     di  = doe(3);
Om = kep_rob(4);     dOm = doe(4);
om = kep_rob(5);     dom = doe(5);
f  = kep_rob(6);      
                     dM  = doe(6);

eta = sqrt(1-e^2);
r = a*(1-e^2)/(1+e*cos(f));
theta = om+f;

% Coordinate relative to the chief orbit: x = x(f)
%                                         y = y(f)
%                                         z = z(f)

x = r/a*da + a*e*sin(f)/eta*dM - a*cos(f)*de;
y = r/eta^3*(1+e*cos(f))^2*dM + r*dom + r*sin(f)/eta^2*(2+e*cos(f))*de + r*cos(i)*dOm;
z = r*(sin(theta)*di - cos(theta)*sin(i)*dOm);
rho = [x y z]';

%--------------------------------------------------------------------------

% Linearized non-dimensional relative orbit motion

% fu = atan(e*dM/(-eta*de));
% fv = fu-pi/2;   % atan(eta*de/(e*dM))
% thetaw = atan(di/(-sin(i)*dOm));
% du = sqrt(e^2*dM^2/eta^2 + de^2);
% dw = sqrt(di^2 + sin(i)^2*dOm^2);

% u = da/a - e*de/(2*eta^2) + du/eta^2*(cos(f-fu)+e/2*cos(2*f-fu));
% v = ((1+e^2/2)*dM/eta^3 + dom + cos(i)*dOm) - du/eta^2*(2*sin(f-fu)+e/2*sin(2*f-fu));
% w = dw*cos(theta-thetaw);
% rho = [u v w]';

%--------------------------------------------------------------------------

return