function [f,I] = dIt_dL(L, L0, P10, P20, th0, I0, dI)

% INPUTS:
% L: longitude of the current state.
% L0: longitude of the initial state.
% P10: initial equinoctial element P1.
% P20: initial equinoctial element P2.
% th0: initial true anomaly
% I0: analytical integrals for the constant tangential case.
%     They are [Ia, Ip1, Ip2] where however Ip1 and Ip2 are not the one given in
%     the Appendix!!!!
% dI:

% OUTPUTS:
% f: integrand of function 31 in the reference
% I:

% Reference: Zuiani, "Extended Analytical Formulas for the Perturbed
% Keplerian Motion Under a Constant Control Acceleration"

% Author: Federico Zuiani
% Comments: Marilena Di Carlo




[a,b]=size(L);

if a>b
    L=L.';
end

% -------------------------------------------------------------------------
% Commented by Federico:
% n=length(L);
% I=zeros(3,n);
% for i=1:n
%     I(:,i)=Integrals_tang_num(L0,L(i),P10,P20);
%
%
% end
% Omom=atan2(P10,P20);
% -------------------------------------------------------------------------

% Since (Omega+omega) is constant:
% L - theta = L0 - theta0
% True anomalies values corresponding to each true longitude value L
th = L - (L0-th0);

% Forward propagation:
if th(end) > th0
    k = (th>pi).*(1+floor((th-pi)/(2*pi)));
    
% Backward propagation
else
    k=(th<-pi).*(-1+ceil((th+pi)/(2*pi)));
end

% -------------------------------------------------------------------------
% Commented by Federico:
% I=Integrals_tang_m(L0,L,P10,P20)+dI*k;
% -------------------------------------------------------------------------


I = I0 + dI*k;
Phi = 1+P10*sin(L)+P20*cos(L);

% Questo dovrebbe essere l'integrando di eq. 31 . Ci sono a, P11 e P21??
% Si. c'era un errore nella formula riportata nel paper. 
f(1,:)=3*(I(1,:)-2*(P10*I(2,:)+P20*I(3,:)))./Phi.^2 - 4 * (1-P10^2-P20^2) * (sin(L).*I(2,:) + cos(L).*I(3,:))./Phi.^3;

return