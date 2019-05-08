% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function K = kepler_K(e, l, P1, P2)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
  This function uses Newton's method to solve Kepler's 
  equation  l = K + P1 cos(K) - P2 sin(K)  
  for the eccentric longitude K,
  given the mean longitude l and P1 and P2

  K  - eccentric longitude (Omega + omega + E) (radians)
  e  - eccentricity, passed from the calling program
  l  - mean longitude (Omega + omega + M)
  P1 - second equinoctial element e * sin(Omega + omega)
  P2 - third equinoctial element e * cos(Omega + omega)
  pi - 3.1415926...

  User m-functions required: none
%}
% ----------------------------------------------

%...Set an error tolerance:
error = 1.e-8;

%...Select a starting value for K:
% if l < pi
%     K = l + e/2;
% else
%     K = l - e/2;
% end
K = l;

%...Iterate until K is determined to within the error tolerance:
ratio = 1;
while abs(ratio) > error
    ratio = (K + P1 * cos(K) - P2 * sin(K) - l)/(1 - P1 * sin(K) - P2 * cos(K));
    K = K - ratio;
end

end %kepler_K
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~