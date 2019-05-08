function  df = dIdrag6_dL3(theta,a,e,incl,omega,coeff, R, J2, Earth_flat)
% Input: theta -> true anomaly
%        a -> semimajor axis
%        e -> eccentricity
%        incl -> inclination
%        omega -> argument of perigee
%        coeff -> coefficients of the chebyshev expansion
%        R -> Earth's radius
%        J2 ->

% Output: df -> integrand of the drag analytic integrals - in this case a
% numerical integration will be realised

% Marilena Di Carlo, 2015
% marilena.di-carlo@strath.ac.uk
% -------------------------------------------------------------------------

n = length(coeff);

p = a * (1-e^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Term of correction for the coupling J2-drag. 
% Reference: Philip Curell, Grace orbit anlaysis tool and parametric
% analysis, pag. 32, eq. 2.15
Delta_r_J2 = J2 * R^2 / (4 * (1-e^2) * a) .* (sin(incl)^2 * cos(2*(omega+theta)) + ...
                                             (3 * sin(incl)^2 - 2) * (1 + e * cos(theta)/(1+sqrt(1-e^2)) + 2 * sqrt(1-e^2) ./ (1 + e * cos(theta)) ) ); 
                                         
h = p ./(1 + e * cos(theta)) - R + Delta_r_J2;

% Remove this if necessary:
fe = 0.00335;
R_flat = R * Earth_flat * (1 - fe * (sin(incl) * sin(omega + theta)).^2) + ...
     R * (1 - Earth_flat);
h = p ./(1 + e * cos(theta)) - R_flat + Delta_r_J2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expansion = 0;

% keyboard
for i = 1 : n
   
    expansion = expansion + coeff(n-i+1,1) * (h).^(i-1);
    
end

% expansion = expansion/(1e-9)*DU^3
% expansion
% keyboard

df = cos(theta) .* sqrt(1+e^2+2*e*cos(theta)) .* expansion ./ ((1+e*cos(theta)).^2);

end