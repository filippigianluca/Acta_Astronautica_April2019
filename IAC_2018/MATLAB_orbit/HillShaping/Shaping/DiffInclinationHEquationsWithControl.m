function [x_prime] = DiffInclinationHEquationsWithControl(u,x,r,G,F_n,flag,i)

%diffEqWithControl: function that writes the differential equations for a
% lt spacecraft with a given control
%
% INPUT
% t: time
% x: [r; v_r; G; u; i; h];
% mu: gravitational parameter
% u: matrix with the control acceleration

% OUTPUT
% xdot: time derivative of cartesian state (r,v). Column vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Column vector

% Write the differential equations
if flag == 0 % Obtain inclination
    
    if length(x) == 2
        x_prime = zeros(2,1) ;
        
        x_prime(1) = r(u) / G(u) * cos(u) * F_n(u) / ( G(u) / r(u)^2 - r(u)*cot(x(1))*sin(u)/G(u) * F_n(u) ) ;
   
        x_prime(2) =  r(u) / G(u) * sin(u) / sin(x(1)) * F_n(u) / ( G(u) / r(u)^2 - r(u)*cot(x(1))*sin(u)/G(u) * F_n(u) ) ;
        
    elseif length(x) == 1
        
        x_prime = r(u) / G(u) * cos(u) * F_n(u) / ( G(u)/r(u)^2 - r(u)*cot(x)*sin(u)/G(u) * F_n(u) ) ;
   
    end

    
elseif flag == 1 && nargin > 6  % Obtain h
    
    x_prime =  r(u) / G(u) * sin(u) / sin(i(u)) * F_n(u) / ( G(u)/r(u)^2 - r(u)*cot(i(u))*sin(u)/G(u) * F_n(u) ) ;

else
    error('Error in Flag or input (not addressed inclination)')
end



end