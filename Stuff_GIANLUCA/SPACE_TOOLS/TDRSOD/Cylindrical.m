%--------------------------------------------------------------------------
% 
% Cylindrical: Computes the fractional illumination of a spacecraft in the 
%              vicinity of the Earth assuming a cylindrical shadow model
% 
% Inputs:
%   r         Spacecraft position vector [m]
%   r_Sun     Sun position vector [m]
%
% Output:
%   nu        Illumination factor:
%             nu=0   Spacecraft in Earth shadow 
%             nu=1   Spacecraft fully illuminated by the Sun
%
% Last modified:   2018/01/04   M. Mahooti
%
%--------------------------------------------------------------------------
function nu = Cylindrical(r, r_Sun)

global const

e_Sun = r_Sun / norm(r_Sun);   % Sun direction unit vector
s     = dot ( r, e_Sun );      % Projection of s/c position 

if ( s>0 || norm(r-s*e_Sun)>const.R_Earth )
    nu = 1.0;
else
    nu = 0.0;
end

