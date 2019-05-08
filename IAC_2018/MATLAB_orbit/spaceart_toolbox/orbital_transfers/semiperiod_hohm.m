function semip = semiperiod_hohm(mu,a)

% semip = semiperiod_hohm(mu,a)
%
% Calculates the semi period of an Hohmann transfer given the semi-major
% axis
%
% INPUT:
%       mu = planetary constant (mu = mass * G) [L^3/T^2]
%       a = semi-major axis
%
% OUTPUT:
%       semip = semi-period of the orbit
%
% - Camilla Colombo - 02/08/2006
%                   - 14/02/2007 - name changed (THohm1.m)
%
% -------------------------------------------------------------------------

semip = pi*sqrt(a^3/mu);
return