function [dI]=Limits_m(P10,P20)

% Function to compute the limits of the analytical integrals of the
% perturbative expansion for the discontinuity at L=pi+2*k*pi, k=1,2,3...
%
% INPUTS;
% P10, P20 : initial values for Equinoctial elements P1 and P2
%
% OUTPUTS:
% dI: 6-row vector containing the values of the "jump" at the
% discontinuity.
% 
% Federico Zuiani

H0=(-1 + P10^2 + P20^2);
H1=sqrt(-H0);
H3=H1^3;
H5=H1^5;
dI11 = (2*pi)/H1;
dI12 = (2*pi)/H3;
dI13 = (pi*(P10^2 + P20^2 + 2))/H5;
dIc2 = -(2*P20*pi)/H3;
dIs2 = -(2*P10*pi)/H3;
dIc3 = -(3*P20*pi)/H5;
dIs3 = -(3*P10*pi)/H5;
dIa2 = P10*pi/(H0*(P20-1));

dI=real([dI11;dI12;dI13;dIc2;dIs2;dIc3;dIs3;dIa2]);

return