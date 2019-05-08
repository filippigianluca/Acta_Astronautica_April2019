% =========================================================================
% Transformation from equinoctial to keplerian elements
% =========================================================================
% Input:  equin  -> equinoctial elements
% Output: kep    -> keplerian elements
% 
% Marilena Di Carlo, March 2014


function kep = eq2kep(equin)


a  = equin(1);
P1 = equin(2);
P2 = equin(3);
Q1 = equin(4);
Q2 = equin(5);
L  = equin(6);   % true longitude


e     = sqrt(P1^2 + P2^2);
incl  = mod(2 * atan(sqrt((Q1^2 + Q2^2))),2*pi);
if incl == 0
    RAAN = 0;
else
    sin_RAAN = Q1 / sqrt(Q1^2 + Q2^2);
    cos_RAAN = Q2 / sqrt(Q1^2 + Q2^2);
    RAAN  = mod(atan2(sin_RAAN, cos_RAAN),2*pi);
end

sin_zeta = P1 / sqrt(P1^2 + P2^2);
cos_zeta = P2 / sqrt(P1^2 + P2^2);
% if abs( P1 / sqrt(P1^2 + P2^2)) > 1
%     keyboard
% end
% if abs(P2 / sqrt(P1^2 + P2^2)) > 1
%     keyboard
% end
% if ~isreal(sin_zeta) || ~isreal(cos_zeta)
%     keyboard
% end
zeta = mod(atan2(sin_zeta, cos_zeta),2*pi);
omega = mod(zeta-RAAN,2*pi);

% omega = mod(atan(P1/P2) - RAAN,2*pi);
true  = mod(L - RAAN - omega,2*pi);

kep = [a, e, incl, RAAN, omega, true];

end