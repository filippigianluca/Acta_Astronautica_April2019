function f = lellipf2(phi,m,errtol)
%
% Exension of lellipf by Thomas Hoffend to calculate the Legendre elliptic
% integral of the first kind for every value of phi
%
% f = lellipf2(phi,m,errtol)
%
% Compute Legendre's (incomplete) elliptic integral F(phi, k).
% Uses a vectorized implementation of Carlson's Duplication Algorithms 
% for symmetric elliptic integrals as found in "Computing Elliptic 
% Integrals by Duplication," by B. C. Carlson, Numer. Math. 33, 1-16 (1979)
% and also found in ACM TOMS Algorithm 577.  Section 4 in the paper cited
% here describes how to convert between the symmetric elliptic integrals
% and Legendre's elliptic integrals.
%
% Returns NaN's for any argument values outside input range.
%
% INPUT:
%       phi = input angle vector size 1xN. There are no constraints on phi
%       m = input parameter vector size 1 or 1xN (m = k^2)
%       errtol = error tolerance for Carlson's algorithms
%
% OUTPUT:
%       f = value of the elliptic integral of the first kind
%
% functions called: lellipf by Thomas Hoffend
%
% evolution of: lellipf
%
% issue - check for which value of m is valid
%
% - Camilla Colombo - 20/02/2007
%
%--------------------------------------------------------------------------

k = sqrt(m);

phi0 = phi;
npi = fix(phi0/(2*pi));
phi = phi0 - npi.*(2*pi); % REM(x,y) is x - n.*y where n = fix(x./y)
% phi = rem(phi0,(2*pi)); % -2*pi<phi<2*pi

% if abs(phi-(phi0-npi*(2*pi)))>=eps(2*pi)
%     fprintf('numerical problem in fix.m and rem.m\nIt seams to be solved - Camilla\n')
%     % phi0
%     % sign(phi0)*eps(phi0)
%     npi = fix((phi0+sign(phi0)*eps(phi0))/(2*pi));
%     phi = rem(phi0,(2*pi));
%     % return
% end

if phi>2*pi
    fprintf('numerical problem in fix.m and rem.m\n - Camilla solve it\n')
end

if (phi>=-pi/2) && (phi<=pi/2)
    f = lellipf(phi,k,errtol);
elseif ((phi>pi/2) && (phi<pi)) || ((phi<-pi/2) && (phi>-pi))
    phisup = -sign(phi)*2*pi+phi;
    n = fix(phisup/(pi));
    esup = 2*n*lellipf(pi/2,k,errtol)+lellipf(phisup-n*pi,k,errtol);
    f = esup+sign(phi)*4*lellipf(pi/2,k,errtol);
elseif (phi>=pi && phi<=3/2*pi) || (phi<=-pi && phi>=-3/2*pi)
    n = fix(phi/(pi));
    f = 2*n*lellipf(pi/2,k,errtol)+lellipf(phi-n*pi,k,errtol);
elseif (phi>3/2*pi) || (phi<-3/2*pi)
    n = fix(phi/(pi));
    phi_pi = phi-n*pi;
    phisup = -sign(phi_pi)*2*pi+phi_pi;
    nsup = fix(phisup/(pi));
    esup = 2*nsup*lellipf(pi/2,k,errtol)+lellipf(phisup-nsup*pi,k,errtol);
    f = 2*n*lellipf(pi/2,k,errtol)+esup+sign(phi_pi)*4*lellipf(pi/2,k,errtol);
% elseif (phi<-pi/2) && (phi>-pi)
%     phisup = 2*pi+phi;
%     n = fix(phisup/(pi));
%     esup = 2*n*lellipf(pi/2,k,errtol)+lellipf(phisup-n*pi,k,errtol);
%     f = esup-4*lellipf(pi/2,k,errtol);
end

f = f + npi*4*lellipf(pi/2,k,errtol);

return    