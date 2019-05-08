function sp = scmb_eq(t,s,u,spl)

% Sun Centered Multi-body equation
%
% sp = scmb_eq(t,s,u,spl)
%
% Equation of the dynamics of the multy-body problem in a Sun-centered
% reference frame, with control
%
% INPUT:
%       t = time [s] since 01/01/2000, 12:00 (MJD2000)
%       s = state vector (position, velocity) [km, km/s]
%       u = control acceleration [km/s^2]
%       spl = vector containing the id numbers of the secondary planets to
%             be considered (< 11):
%               1:	Mercury
%             	2:	Venus
%            	3:	Earth
%             	4:	Mars
%             	5:	Jupiter
%             	6:	Saturn
%             	7:	Uranus
%            	8:	Neptune
%            	9:  Pluto
%            	10: Moon
%
% OUTPUT:
%       sp = derivative of the state
%
% FUNCTIONS CALLED: astro_constants, EphSS
%
% - Camilla Colombo - 11/05/2007
%  Revised by Matteo Ceriotti - 13-05-2007
%
% ------------------------- - SpaceART Toolbox - --------------------------

% Sun centered vector
r = s(1:3);
v = s(4:6);
u = u(:);

% Primary body
mu_sun = astro_constants(4);

% Secondaries bodies (it can not be the Sun, if you want the Sun another
% function in place of uplanet should be used)

a = zeros(3,1);
for i = 1:length(spl)
    mu = astro_constants(spl(i)+10);
    rp = EphSS(spl(i),t/86400);
    if i == 10
        rpE = EphSS(3,t/86400);
        rp = rp'+rpE';
    end
    rpl = r - rp'; % position vector from planet i
    nrpl = sqrt(rpl(1)^2+rpl(2)^2+rpl(3)^2);
    a = a - mu(i)/nrpl^3*rpl;
end

sp = [       v        ;...
      -mu_sun/sqrt(r(1)^2+r(2)^2+r(3)^2)^3*r + a + u];

return