function f = lellipf3_mod(phi,m,errtol)
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
% Modified by Marilena Di Carlo - 16/04/2015
% Scrivere dove
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

f=0*phi;
logi_a=(phi>=-pi/2) & (phi<=pi/2);
logi_b=((phi>pi/2) & (phi<pi)) | ((phi<-pi/2) & (phi>-pi));
logi_c=(phi>=pi & phi<=3/2*pi) | (phi<=-pi & phi>=-3/2*pi);
logi_d=~(logi_a|logi_b|logi_c);

if any(logi_a)
    f(logi_a) = lellipf(phi(logi_a),k(logi_a),errtol);
end
if any(logi_b|logi_c|logi_d)
    Df=lellipf(pi/2,k,errtol);
    if any(logi_b)
        phisup = -sign(phi(logi_b))*2*pi+phi(logi_b);
        n = fix(phisup/(pi));
        esup = 2.*n.*Df(logi_b)+lellipf(phisup-n*pi,k(logi_b),errtol);
        f(logi_b) = esup+sign(phi(logi_b)).*4.*Df(logi_b);
    end
    if any(logi_c)
        n = fix(phi(logi_c)/(pi));
        f(logi_c) = 2.*n.*Df(logi_c)+lellipf(phi(logi_c)-n*pi,k(logi_c),errtol);
    end
    if any(logi_d)
        n = fix(phi(logi_d)/(pi));
        phi_pi = phi(logi_d)-n*pi;
        phisup = -sign(phi_pi)*2*pi+phi_pi;
        nsup = fix(phisup/(pi));
        esup = 2*nsup*Df+lellipf(phisup-nsup*pi,k,errtol);
        f(logi_d) = 2*n*Df+esup+sign(phi_pi)*4*Df;
    end
elseif npi
    Df=lellipf(pi/2,k,errtol);
else 
    Df=0;
end

f = f + npi.*4.*Df;

return    