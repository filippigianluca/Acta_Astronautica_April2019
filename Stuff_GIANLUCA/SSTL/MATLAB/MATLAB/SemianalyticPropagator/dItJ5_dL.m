function [f,I] = dItJ5_dL(L, L0, a0, P10, P20, Q10, Q20, th0, a1_J5, P11_J5, P21_J5)




[a,b]=size(L);

if a>b
    L=L.';
end


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



% I = I0 + dI*k;
Phi = 1+P10*sin(L)+P20*cos(L);


B0 = sqrt(1 - P10^2 - P20^2);
G0 = 1 + Q10^2 + Q20^2;

% a1 = a0 / B0^12 / G0 * I_J5_a;
% 
% P11 = - I_J5_P11 / (2 * G0 * B0^10) + P20 / (B0^10 * G0) * (1 - Q10^2 - Q20^2) * I_J5_P12;
%  
% P21 = - I_J5_P21 / (2 * G0 * B0^10) - P10 / (8 * B0^10 * G0) * (1 - Q10^2 - Q20^2) * I_J5_P12;

f(1,:)= 3/2 * B0^2 * a1_J5 ./ Phi.^2 - a0 * ( (3 * P10 ./ Phi.^2 + 2 * B0^2 * sin(L) ./ Phi.^3) * P11_J5 + ...
                                              (3 * P20 ./ Phi.^2 + 2 * B0^2 * cos(L) ./ Phi.^3) * P21_J5);

return