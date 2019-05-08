function dv=dv_insertion(s,t,vf,rp,e)
% Computes the dv to insert in an orbit around a planet.
%   
%   dv = dv_insertion(s, t, vf, rp, e)
%
% It assumes a patched conics framework, and performs the manoeuvre at the
% pericentre of the incoming hyperbola to achieve the target orbit.
%
% INPUT
%   s   Planet (1 <= s <= 9).
%   t   Time of the insertion (final time of the trajectory) [d, MJD2000].
%   vf  Absolute cartesian velocity at the planet arrival [km/s].
%   rp  Radius of the pericentre of the orbit to achieve around the planet
%       [km].
%   e   Eccentricity of the orbit to achieve.
%
% OUTPUT
%   dv  Delta-v needed to insert into the defined orbit [km/s].
%
% FUNCTIONS CALLED: astro_constants, EphSS
%
% Matteo Ceriotti, 19-05-2007
% Revised by Camilla Colombo, 21-05-2007
%
% ------------------------- - SpaceART Toolbox - --------------------------

vf=vf(:);

mu=astro_constants(10+s);
[xp2,vp2]=EphSS(s,t);
vp2=vp2';

vinf=norm(vf-vp2); % Velocity relative to the planet (at infinite on the hyperbola)
vp1=sqrt(vinf^2+2*mu/rp); % Velocity at the pericentre of the hyperbola
vp2=sqrt(2*mu/rp-mu/rp*(1-e)); % Velocity at the pericentre of the target orbit
dv=abs(vp1-vp2); % Delta-v at the pericentre

return