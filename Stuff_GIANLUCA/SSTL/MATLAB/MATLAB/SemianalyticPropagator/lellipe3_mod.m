function e = lellipe3_mod(phi,m,errtol)
%
% Exension of lellipe by Thomas Hoffend to calculate the Legendre elliptic
% integral of the second kind for every value of phi
%
% e = lellipe2(phi,m,errtol)
%
% Compute Legendre's (incomplete) elliptic integral E(phi,k).
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
%       e = value of the elliptic integral of the second kind
%
% functions called: lellipe by Thomas Hoffend
%
% evolution of: lellipe
%
% issue - check for which value of m is valid
%
% - Camilla Colombo - 20/02/2007
% Modified by Marilena Di Carlo - 16/04/2015
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

e=0*phi;
logi_a=(phi>=-pi/2) & (phi<=pi/2);
logi_b=((phi>pi/2) & (phi<pi)) | ((phi<-pi/2) & (phi>-pi));
logi_c=(phi>=pi & phi<=3/2*pi) | (phi<=-pi & phi>=-3/2*pi);
logi_d=~(logi_a|logi_b|logi_c);

if any(logi_a)
    e(logi_a) = lellipe(phi(logi_a),k(logi_a),errtol);
end
if any(logi_b|logi_c|logi_d)
    De=lellipe(pi/2,k,errtol);
    if any(logi_b)
        phisup = -sign(phi(logi_b))*2*pi+phi(logi_b);
        n = fix(phisup/(pi));
        esup = 2.*n.*De(logi_b)+lellipe(phisup-n*pi,k(logi_b),errtol);
        e(logi_b) = esup+sign(phi(logi_b)).*4.*lellipe(pi/2,k(logi_b),errtol);
    end
    if any(logi_c)
        n = fix(phi(logi_c)/(pi));
        e(logi_c) = 2.*n.*De(logi_c)+lellipe(phi(logi_c)-n*pi,k(logi_c),errtol);
    end
    if any(logi_d)
        n = fix(phi(logi_d)/(pi));
        phi_pi = phi(logi_d)-n*pi;
        phisup = -sign(phi_pi)*2*pi+phi_pi;
        nsup = fix(phisup/(pi));
        esup = 2*nsup*De+lellipe(phisup-nsup*pi,k,errtol);
        e(logi_d) = 2*n*De+esup+sign(phi_pi)*4*De;
    end
elseif npi
    De=lellipe(pi/2,k,errtol);
else
    De=0;
end

e = e + npi.*4.*De;

return    