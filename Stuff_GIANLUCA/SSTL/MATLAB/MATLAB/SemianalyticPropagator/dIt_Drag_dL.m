function [f,I] = dIt_Drag_dL(L, L0, P10, P20, th0, IDrag, dI)

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

% Author:  Marilena Di Carlo




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

Omom = atan2(P10,P20);


% -------------------------------------------------------------------------
% Secondo me questa parte non serve perche' limiti destro e sinistro degli
% integrali nella discontinuita` in pi sono uguali. Attendo conferma da
% Federico e poi lo cancello
% Forward propagation:
if th(end) > th0
    k = (th>pi).*(1+floor((th-pi)/(2*pi)));
    
% Backward propagation
else
    k=(th<-pi).*(-1+ceil((th+pi)/(2*pi)));
end

% I = IDrag + dI*k;
% -------------------------------------------------------------------------

Phi = 1+P10*sin(L)+P20*cos(L);
e0 = sqrt(P10^2 + P20^2);
B = sqrt(1 - P10^2 - P20^2);

P10_vec = P10 * ones(1,length(L));
P20_vec = P20 * ones(1,length(L));

f = 1.5 * B^2 * (e0^2 * IDrag(1,:) + IDrag(2,:)) ./ Phi.^2 + ...
   - (B^2)/2 * (3 * P10_vec / Phi.^2 + 2 * B^2 * sin(L) ./ Phi.^3 ) .* (sin(Omom) * (e0 * IDrag(1,:) + e0 * IDrag(4,:) + 2 * IDrag(6,:) + e0 * IDrag(7,:)) + 2 * cos(Omom) * IDrag(5,:) )  + ...
   - (B^2)/2 * (3 * P20_vec / Phi.^2 + 2 * B^2 * cos(L) ./ Phi.^3 ) .* (cos(Omom) * (e0 * IDrag(1,:) + e0 * IDrag(4,:) + 2 * IDrag(6,:) + e0 * IDrag(7,:)) - 2 * sin(Omom) * IDrag(5,:) ) ;

end