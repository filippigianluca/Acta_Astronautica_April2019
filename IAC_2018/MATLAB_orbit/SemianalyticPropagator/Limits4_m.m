function [dI]=Limits4_m(P10,P20)

% Function to compute the limits of the analytical integrals of the
% perturbative expansion for the discontinuity at L=pi+2*k*pi, k=1,2,3...
%
% INPUTS;
% P10, P20 : initial values for Equinoctial elements P1 and P2
%
% OUTPUTS:
% dI: 6-row vector containing the values of the "jump" at the
% discontinuity.

% Federico Zuiani, Marilena Di Carlo
% marilena.di-carlo@strath.ac.uk

H0 = (-1 + P10^2 + P20^2);
H1 = sqrt(-H0);
H3 = H1^3;
H5 = H1^5;
 
dI11 = (2*pi)/H1;
dI12 = (2*pi)/H3;
dI13 = (pi*(P10^2 + P20^2 + 2))/H5;

dIc2 = -(2*P20*pi)/H3;
dIs2 = -(2*P10*pi)/H3;

dIc3 = -(3*P20*pi)/H5;
dIs3 = -(3*P10*pi)/H5;

% Chi e' questo?
dIa2 = P10*pi/(H0*(P20-1));

dIcc3 = ((1 - P10^2 + 2*P20^2)*pi)/H5;
dIcs3 = (3*P10*P20*pi)/H5;

dIss3 = ((1 + 2*P10^2 - P20^2)*pi)/H5;

% Controllare questi: ma dove vengono usati?
dI3c3=(P20*pi*(-2*(3*P10^2 - P20^2) + (9*P10^6 - 2*P20^2 + 5*P20^4 - 6*P20^6 + 3*P10^4*(-5 + 4*P20^2) + P10^2*(6 - 10*P20^2 - 3*P20^4))/H5))/(P10^2 + P20^2)^3;
dI2c1s3=(P10*pi*(-2*(P10^2 - 3*P20^2) + (3*P10^6 - P10^4*(5 + 6*P20^2) + P10^2*(2 + 10*P20^2 - 21*P20^4) - 3*P20^2*(2 - 5*P20^2 + 4*P20^4))/H5))/(P10^2 + P20^2)^3;
dI1c2s3=(P20*pi*(2*(3*P10^2 - P20^2) - (12*P10^6 - 2*P20^2 + 5*P20^4 - 3*P20^6 + 3*P10^4*(-5 + 7*P20^2) + 2*P10^2*(3 - 5*P20^2 + 3*P20^4))/H5))/(P10^2 + P20^2)^3;
dI3s3=(P10*pi*(2*(P10^2 - 3*P20^2) - (6*P10^6 - 6*P20^2 + 15*P20^4 - 9*P20^6 + P10^4*(-5 + 3*P20^2) + 2*P10^2*(1 + 5*P20^2 - 6*P20^4))/H5))/(P10^2 + P20^2)^3;
dI4c3=-((3*pi*(2*(P10^4 - 6*P10^2*P20^2 + P20^4) + (P10^4*(-2 + P10^2)*(-1 + P10^2)^2 - P10^2*(-12 + P10^2)*(-1 + P10^2)^2*P20^2 - (2 + 25*P10^2 - 36*P10^4 + 9*P10^6)*P20^4 + (5 + 14*P10^2 - 11*P10^4)*P20^6 - 4*(1 + P10^2)*P20^8)/H5))/(P10^2 + P20^2)^4);
dI3c1s3=(3*P10*P20*pi*(8*(P10 - P20)*(P10 + P20) + (3*P10^8 + 8*P20^2 - 20*P20^4 + 15*P20^6 - 2*P20^8 + P10^6*(-15 + 7*P20^2) + P10^2*(-8 + 15*P20^4 - 3*P20^6) + P10^4*(20 + 3*P20^2*(-5 + P20^2)))/H5)/(P10^2 + P20^2)^4);
dI2c2s3=(pi*(6*(P10^4 - 6*P10^2*P20^2 + P20^4) + (2*P10^10 - P10^8*(11 + 5*P20^2) + P10^6*(15 + 46*P20^2 - 25*P20^4) + (-1 + P20)*P20^4*(1 + P20)*(6 - 9*P20^2 + 2*P20^4) + P10^2*P20^2*(36 - 75*P20^2 + 46*P20^4 - 5*P20^6) - P10^4*(6 + 75*P20^2 - 114*P20^4 + 25*P20^6))/H5)/(P10^2 + P20^2)^4);
dI1c3s3=((3*P10*P20*pi*(-8*(P10 - P20)*(P10 + P20) - (2*P10^8 + 8*P20^2 - 20*P20^4 + 15*P20^6 - 3*P20^8 + 3*P10^6*(-5 + P20^2) + P10^2*(-8 + 15*P20^4 - 7*P20^6) + P10^4*(20 - 3*P20^2*(5 + P20^2)))/H5))/(P10^2 + P20^2)^4);
dI4s3=-((3*pi*(2*(P10^4 - 6*P10^2*P20^2 + P20^4) - (P10^2*P20^2*(-12 + P20^2)*(-1 + P20^2)^2 - P20^4*(-2 + P20^2)*(-1 + P20^2)^2 + 4*P10^8*(1 + P20^2) + P10^6*(-5 - 14*P20^2 + 11*P20^4) + P10^4*(2 + 25*P20^2 - 36*P20^4 + 9*P20^6))/H5))/(P10^2 + P20^2)^4);


% In the file calling this function three additional zero elements will be
% placed before the vector dI, in such a way as that dI11 will be the 24th
% elements, as it should be (refer to the Guide to the Analytic Integrals)
dI=[zeros(20,1);[dI11; dI12; dI13; dIc2; dIs2; dIc3; dIs3; dIa2; ...
                 dIcc3; dIcs3; dIss3; dI3c3; dI2c1s3; dI1c2s3; dI3s3; ...
                 dI4c3; dI3c1s3; dI2c2s3; dI1c3s3; dI4s3]];

return