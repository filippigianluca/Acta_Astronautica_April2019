function ts=localsolartime(X,E,Alpha,t,one_day)
% ts=(t/0.00274-floor(t/0.00274))*2*pi;
% ts=(((t / 0.00274) - floor(t / 0.00274)) * 86400);

% orbital elements
a=X(1);             % semimajor axis
e=X(2);             % eccentricity
I=X(3);             % inclination
w=X(4);             % argurment of perigee
CapitalOmega=X(5);  % RAAN
f=X(6);             % true anomaly

% ts=asin(((-1)+e.*cos(E)).^(-1).*(1+(-1).*((-1)+e.*cos(E)).^(-2).*sin(I) ...
%   .^2.*((1+(-1).*e.^2).^(1/2).*cos(w).*sin(E)+((-1).*e+cos(E)).*sin(w)) ...
%   .^2).^(-1/2).*(cos(I).*cos(Alpha+(-1).*CapitalOmega).*((1+(-1).*e.^2).^(1/2).*cos(w).* ...
%   sin(E)+((-1).*e+cos(E)).*sin(w))+((e+(-1).*cos(E)).*cos(w)+(1+(-1).* ...
%   e.^2).^(1/2).*sin(E).*sin(w)).*sin(Alpha+(-1).*CapitalOmega)));

ts=(t/one_day-floor(t/one_day))*24;

end
