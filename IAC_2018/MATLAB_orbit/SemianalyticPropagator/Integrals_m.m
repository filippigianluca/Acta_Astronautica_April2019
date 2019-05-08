function [I,I_test]=Integrals_m(L0, L, P10, P20)

% Calculates analytical integrals for the perturbative expansion.
%
% INPUTS:
% L0: longitude of the initial state.
% L: longitude of the current state.
% P10, P20: values at x(1) for Equinoctial elements P1 and P2.
%
% OUTPUTS:
% I: 6-row vector with the integrals between x(1) and x(2).

% Reference: Zuiani, "Extended Analytical Formulas for the Perturbed
% Keplerian Motion Under a Constant Control Acceleration"

% Author: Federico Zuiani
% Comments and code tidying: Marilena Di Carlo



sL  = sin(L);
sL0 = sin(L0);
cL  = cos(L);
cL0 = cos(L0);


% First condition: Omega + omega is NOT approximately 90 
if ~(abs(1/P20)>1e16)
    
    % =====================================================================
    % Equation [3]: B = sqrt(1 - P10^2 - P20^2)
    % H1 is B
    H0 = (-1 + P10^2 + P20^2);
    H1 = sqrt(-H0);
    H3 = H1^3;
    H5 = H1^5;

    
    
    % =====================================================================
    % Equation [3]: Phi = 1 + P10 sinL + P20 cosL
    M1 = (1 + P20*cL + P10*sL);
    M2 = M1^2;
    M10 = (1 + P20*cL0 + P10*sL0);
    M20 = M10^2;

    
    
    % =====================================================================
    % Equation [10] in Appendix (Lambda0)
    K = atan( (-P10 + (-1 + P20)*tan(L/2)) / (H1) );
    K0 = atan( (-P10 + (-1 + P20)*tan(L0/2)) / (H1) );
    K1=(K - K0);

    
    % ---------------------------------------------------------------------
    % Commented by Federico:
%     S=(P10 + sL);
%     R=(P10*(3*M1-H0)+sL*(3*M1+2*M1*H0-H0));
%     S0=(P10 + sL0);
%     R0=(P10*(3*M10-H0)+sL0*(3*M10+2*M10*H0-H0));
    % ---------------------------------------------------------------------
    
    % =====================================================================
    % Equation [1] in Appendix
    I11 = -(2 * K1 / H1);

    
    
    % =====================================================================
    % Equation [2] in Appendix 
    % Secondo me ha i segni sbagliati. Non dovrebbe essere:
    % I12 = -2/H3 * K1 - 1/(P20*H0)*( (P10 + (P10^2 + P20^2) * sL) / (M1) - ( P10 + (P10^2 + P20^2)*sL0) / (M10));
    % No, perche' usa H0, che e' l'inverso di B. 
    I12 = -2/H3 * K1 + 1/(P20*H0)*( (P10 + (P10^2 + P20^2) * sL) / (M1) - ( P10 + (P10^2 + P20^2)*sL0) / (M10));

    
    
    
    % =====================================================================
    % Equation [3] in Appendix
    % Anche qui problemi di segni?!?!?
%     I13 = (1/2) * (-2*(2 + P10^2 + P20^2)/H5*K1...
%         - 1/(P20*H0) * (( (P10 + (P10^2 + P20^2)*sL)/M2 - (P10*(2 + P10^2 + P20^2) + 3*(P10^2 + P20^2)*sL)/(H0*M1)) - ((P10 + (P10^2 + P20^2)*sL0)/M20 - (P10*(2 + P10^2 + P20^2) + 3*(P10^2 + P20^2)*sL0)/(H0*M10))));
    % No, perche' usa H0, che e' l'inverso di B. 
    I13 = (1/2) * (-2*(2 + P10^2 + P20^2)/H5*K1...
        + 1/(P20*H0) * (( (P10 + (P10^2 + P20^2)*sL)/M2 - (P10*(2 + P10^2 + P20^2) + 3*(P10^2 + P20^2)*sL)/(H0*M1)) - ((P10 + (P10^2 + P20^2)*sL0)/M20 - (P10*(2 + P10^2 + P20^2) + 3*(P10^2 + P20^2)*sL0)/(H0*M10))));
    % ---------------------------------------------------------------------
    % Commented by Federico:
    %         + 1/P20/H0^2*((-R+sL*H0*(H0-M1))/M2-(-R0+sL0*H0*(H0-M10))/M20)-P10/P20/H0^2);
    % ---------------------------------------------------------------------
    
    
    
    % =====================================================================
    % Equation [4] in Appendix
    Ic2 = (2*P20*K1/H3) - ( (P10 + sL)/M1 - (P10 + sL0)/M10) / H0;
    % ---------------------------------------------------------------------
    % Commented by Federico:
    %     (S/M1-S0/M10)/H0;
    % ---------------------------------------------------------------------

    
    
    
    % =====================================================================   
    % Equation [5] in Appendix
    Is2 = (2*P10*K1/H3) + ( (-1 + P20^2 - P10*sL)/M1 - (-1 + P20^2 - P10*sL0)/M10 ) / (P20*H0);
    % ---------------------------------------------------------------------
    % Commented by Federico:
    %     ((1 - P10*S/H0)/M1 - (1 - P10*S0/H0)/M10)/P20;
    % ---------------------------------------------------------------------

    
    
    
    % =====================================================================   
    % Equation [6] in Appendix
    Ic3 = (1/2)*((6*P20/H5*K1)...
        + 1/H0 * ( ( - (P10 + sL)/(M2) + (3*P10 + (1 + 2*P10^2 + 2*P20^2)*sL)/(H0*M1) ) -...
        (- (P10 + sL0)/(M20) + (3*P10 + (1 + 2*P10^2 + 2*P20^2)*sL0)/(H0*M10))));
    % ---------------------------------------------------------------------
    % Commented by Federico:
    %     + (R/M2-R0/M20)/H0^2);
    % ---------------------------------------------------------------------

    
    
    
    
    % =====================================================================   
    % Equation [7] in Appendix
    Is3 = (1/2)*((6*P10/H5*K1)...
        + 1/(P20*H0)*(((-1 + P20^2 - P10*sL)/M2 + (P10*(3*P10 + (1 + 2*P10^2 + 2*P20^2)*sL))/(H0*M1)) -...
        ((-1 + P20^2 - P10*sL0)/M20 + (P10*(3*P10 + (1 + 2*P10^2 + 2*P20^2)*sL0))/(H0*M10))));
    % ---------------------------------------------------------------------
    % Commented by Federico:
    %         + P10/P20*((R/H0^2+1)/M2-(R0/H0^2+1)/M20));
    % ---------------------------------------------------------------------
    
    
    
    
    % ===================================================================== 
    % Cos'e'?
    Ia2=-(1/(2*(-1 + P20)*H3))*(...
        2*(-1 + P20)*(K^2-K0^2) -...
        (P10^2 + (-1 + P20)*P20)*(cos(2*K)-cos(2*K0)) -...
        P10*H1*(sin(2*K)-sin(2*K0)) +...
        2*(P10*H1*(K*cos(2*K)-K0*cos(2*K0)) - (P10^2 + (-1 + P20)*P20)*(K*sin(2*K)-K0*sin(2*K0))));
    

    
    %%% CANCELLARE
%     I_test = I11 * I12 - 2 /(H1^4) *(K^2-K0^2) + 1/H1^2 *(1/M1 - 1/M10) + 2 * P10 / (P20 * H1^3)* (K-K0);
    
    
    
% Second condition: (Omega + omega) ~= 0 
% If the first condition does not applies it means that
% (Omega + omega) = 90 but still e~=0. 
% This means that P2 = e cos(Omega + omega) = 0.
% In the following equations therefore P2 will be removed. 
elseif ~(abs(1/P10)>1e16)
    
    H1 = sqrt(1 - P10^2);
    H3 = H1^3;
    H5 = H1^5;
    
    % Non dovrebbe essere con il segno meno? Si, ma guardando alla
    % definizione di I11 #dove non c'e' il meno) tutto torna!
    K =  atan((P10 + tan(L/2))/H1);
    K0 = atan((P10 + tan(L0/2))/H1);
    K1=(K - K0);
    
    M1 = (1 + P10*sL);
    M2 = M1^2;
    M10 = (1 + P10*sL0);
    M20 = M10^2;
    
    I11 = 2*K1 / H1;
    I12 = (2*K1 + (P10*H1 * (cL/M1 - cL0/M10)))/H3;
    I13 = 1/(2*(1 - P10^2)^2) * ( (2*(2 + P10^2)*K1)/H1 + P10 * (cL * (4 - P10^2 + 3*P10*sL)/M2 - cL0*(4 - P10^2 + 3*P10*sL0)/M20) );
    
    Ic2 = -(1/M1-1/M10)/P10;
    Is2 = -((cL/M1 - cL0/M10)/(1 - P10^2) + 2*P10*K1/H3);
    
    Ic3 = 1/(2*P10)*(-1/M2+1/M20);
    
    % Il seguente integrale non mi torna:
    Is3 = 1/(2*H5)*(-6*P10*K1 - H1*(cL*(2 + P10^2 + (P10 + 2*P10^3)*sL)/M2 - cL0*(2 + P10^2 + (P10 + 2*P10^3)*sL0)/M20));
    
    % Cos'e Ia2?
    Ia2=-(2*(K^2-K0^2) + 2*P10*(H1*(K*cos(2*K)-K0*cos(2*K0)) + P10*(K*sin(2*K)-K0*sin(2*K0))) + P10*(P10*(cos(2*K)-cos(2*K0)) - H1*(sin(2*K)-sin(2*K0))))/(2*H3);
    
%     %%% CANCELLARE
%     I_test = I11 * I12 - 2 /(H1^4) *(K^2-K0^2) + 1/H1^2 *(1/M1 - 1/M10) + ...
%         2 * P10 / (P20 * H1^3)* (K-K0);
    
    
% Third condition: P1 and P2 are both zero
else
    
    % In this case, Phi = 1 + P1 sinL + P2 cosL = 1
    % I11, I12 and I13 have integrands equal to 1
    I11 = L-L0;
    I12 = I11;
    I13 = I11;
    
    % The integrand of Ic2 is cos L
    Ic2 = sL-sL0;
    
    % The integrand of Is2 is sin L
    Is2 = cL0-cL;
    
    Ic3 = Ic2;
    Is3 = Is2;
    
    
    % Cos'e Ia2?
    LL = mod(L+pi,2*pi)-pi;
    LL0 = mod(L0+pi,2*pi)-pi;
    Ia2 = (LL^2-LL0^2)/4 - (LL*atan(tan(LL/2))-LL0*atan(tan(LL0/2)));
end



% The output of this file are elements 24 to 31 in the bigger vector of
% analytic integrals (refer to AnalyticIntegralsGuide.m)
I = [I11; I12; I13; Ic2; Is2; Ic3; Is3; Ia2];

% keyboard

return