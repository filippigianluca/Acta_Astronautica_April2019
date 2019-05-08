function [e,w] = Geop_effect_low_e(X, constant,aph, A, t, to, one_day)
% Since Kuala analysis is unsuitable for low eccentric and inclination 
% values and we use nonsingular elements h and k (equinoctical elements)

% orbital elements
a=X(1);             % semimajor axis
e=X(2);             % eccentricity
I=X(3);             % inclination
w=X(4);             % argurment of perigee
CapitalOmega=X(5);  % RAAN
f=X(6);             % true anomaly

% propagation time
tp=propagation_time(t,to,one_day);

% aph = -1.457396617476020;   % derived from the setting time to zeros
% A = -0.002455327748045;     % derived from the setting time to zeros

B = -0.5 * (constant.Rearth/a) * (constant.J3/constant.J2) * sin(I);
% B=(1/512).*a.^(-7).*constant.J2.^(-1).*constant.Rearth.*(a.^(-3).*constant.mu).^(-1/2).*(3+5.*cos(2.*I)) ... % Derived from Cook 1966 papaer
%   .^(-1).*sin(I).*((-256).*a.^6.*constant.J3.*(3+5.*cos(2.*I))+80.*a.^4.*constant.J5.* ...
%   constant.Rearth.^2.*(15+28.*cos(2.*I)+21.*cos(4.*I))+140.*a.^2.*constant.J7.*constant.Rearth.^4.*((-64)+ ...
%   432.*sin(I).^2+(-792).*sin(I).^4+429.*sin(I).^6)+105.*constant.J9.*constant.Rearth.^6.*( ...
%   128+11.*sin(I).^2.*((-128)+416.*sin(I).^2+(-520).*sin(I).^4+221.* ...
%   sin(I).^6)));




sigma = 3 * sqrt(constant.mu/a^3) * constant.J2 * (constant.Rearth/a)^2 * (1 - 5/4 * sin(I)^2);

k = A * cos(sigma*tp + aph);

h = A * sin(sigma*tp + aph) + B;



e=sqrt((h^2) + (k^2));
w=atan2(h, k);

end
