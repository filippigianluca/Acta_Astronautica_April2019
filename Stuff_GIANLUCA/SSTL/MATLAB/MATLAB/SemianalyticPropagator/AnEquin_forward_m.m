function [Equin,t]=AnEquin_forward_m(L,Equin0,epsilon,alfa,beta,muadim)

% Function which implements the analytical formulas for the perturbative
% expansion of Keplerian motion (no perturbations, except for the low thrust
% one, but possible to realise backward propagation)

% INPUTS:
% L: row vector containing the inital value for longitude L0 and the points
%   in which the Equinoctial elements should be evaluated.
% Equin0: column vector containing the initial state (at L0) in Equinoctial
%   elements.
% epsilon,alfa,beta: modulus, azimuth and elevation (in radians) of the
%   perturbing acceleration vector, in the radial, transverse, normal
%   reference system.
% muadim: gravitational constant of the central body, with dimensions
%   consistent with the other quantities.
%
% OUTPUTS:
% Equin: Matrix containing the column vectors in Equinoctial elements for
%   each value of longitude contained in L.
%
% Calls "Integrals." and "Limits.m"
%
%   Federico Zuiani 2010
% Overloaded version

%% Initialization
if (L(1)==Equin0(6))&&(length(L)>1)
    L=L(2:end);
end
n=length(L);
L0=Equin0(6);
a0=Equin0(1);
P10=Equin0(2);
P20=Equin0(3);
Q10=Equin0(4);
Q20=Equin0(5);
% l0=Equin0(6);
e0=sqrt(P10^2+P20^2);
B=(1-P10^2-P20^2)^(1/2);
p0=a0*B^2;
h0=sqrt(muadim*p0);
h0=sqrt(1/muadim*p0);
if e0>1
%     warning('Cannot treat open trajectories')
end

% L0 should be comprised between [0,2*pi)
if L(end)>L0
    while L0>=pi
        L0=L0-2*pi;
    end
    while L0<-pi
        L0=L0+2*pi;
    end
    Ic=Integrals(L0,pi,P10,P20);
else
    while L0>pi
        L0=L0-2*pi;
    end
    while L0<=-pi
        L0=L0+2*pi;
    end
    Ic=Integrals(L0,-pi,P10,P20);
end

% Adjust the rest of L accordingly
L=L-(Equin0(6)-L0);

% Calculation of limits for the discontinuity at L=pi+2*k*pi, k=1,2,3...
dI=Limits_m(P10,P20);


%% Begin cycle

Ibas = zeros(8,1);
corr=0;
Equin=NaN*ones(6,n);
t=NaN*ones(1,n);
% Adjust the initial value for multiplier n_corr which accounts for the
% discontinuity.
n_corr=1;
% Equin(:,1)=[a0 P10 P20 Q10 Q20 L0]';
for i=1:n
    
    % Discontinuity correction
    if e0~=0 % Correction is needed only for elliptical orbits
        if L(end)>L0 % Case of increasing L (forward propagation)
            while L(i)>(n_corr*pi)
                Ibas=Ibas+dI;
                n_corr=n_corr+2;
                corr=corr+dI(1)*Ic(2);
                Ic(2)=Ic(2)+dI(2);
            end
        else % Case of decreasing L (backward propagation)
            while L(i)<-(n_corr*pi)
                Ibas=Ibas-dI;
                n_corr=n_corr+2;
                corr=corr-dI(1)*Ic(2);
                Ic(2)=Ic(2)-dI(2);
            end
        end
    end
    % Evaluate analytical integrals
    I=Ibas+Integrals_m(L0,L(i),P10,P20);
    
    % Calculate first order variation of Equinoctial elements
    a11 = 2*(h0*a0/muadim)^2*cos(beta)*(cos(alfa)*(P20*I(5)-P10*I(4))+sin(alfa)*I(1));
    P11=h0^4/muadim^3*(cos(beta)*(-cos(alfa)*I(4)+sin(alfa)*(P10*I(3)+I(7)+I(5)))+sin(beta)*...
        P20*(-Q10*I(6)+Q20*I(7)));
    P21=h0^4/muadim^3*(cos(beta)*(cos(alfa)*I(5)+sin(alfa)*(P20*I(3)+I(6)+I(4)))+sin(beta)*...
        P10*(Q10*I(6)-Q20*I(7)));
    Q11=h0^4/(2*muadim^3)*(1+Q10^2+Q20^2)*sin(beta)*I(7);
    Q21=h0^4/(2*muadim^3)*(1+Q10^2+Q20^2)*sin(beta)*I(6);
    if nargout>1
        t00=h0^3*I(2);
        % Standard
%         t11s=((2*a0*(2*P10^2*P21-3*P10*P11*P20-P21*(P20^2+2))-3*a11*P20*(P10^2+P20^2-1))*I(6)/2-(2*a0*(P10^2*P11+3*P10*P20*P21+2*P11*(1-P20^2))+3*a11*P10*(P10^2+P20^2-1))*I(7)/2-3*(2*a0*(P10*P11+P20*P21)+a11*(P10^2+P20^2-1))*I(3)/2)*sqrt(a0*(1-P10^2-P20^2));
%         t11s=3/2*sqrt(a0/muadim)*B^3*I(2)*a11 -...
%             sqrt(a0^3/muadim)*B*(2*B^2*I(7)+3*P10*I(2))*P11 -...
%             sqrt(a0^3/muadim)*B*(2*B^2*I(6)+3*P20*I(2))*P21;
%         % Numerical
%         Ixx=Integrals_t_num(L0,L(i),P10,P20,Q10,Q20,alfa,beta);
%         t11n=3*a0^(7/2)*B^5/(muadim^(3/2))*cos(beta)*(cos(alfa)*Ixx(1)+sin(alfa)*Ixx(2))-...
%             a0^(7/2)*B^5/(muadim^(3/2))*(2*B^2*Ixx(3)+3*Ixx(4));
        
        % Modified analytic
%         I0=Integrals01(L0,P10,P20);
        if e0~=0
%             keyboard
            t11=3*a0^(7/2)*B^5/(muadim^(3/2)) * cos(beta) * (cos(alfa)*(I(3)-I(2)/(1+P10*sin(L0)+P20*cos(L0)))...
                -2*sin(alfa)/B*(I(8)-(I(2)*atan((-P10+(P20-1)*tan(L0/2))/B)+B/2*(Ibas(1)*I(2)-corr))));
            
            
%             keyboard
            % CANCELLARE
%             [~,I_test] = Integrals_m(L0,L(i),P10,P20);
%               3*a0^(7/2)*B^5/(muadim^(3/2)) * cos(beta) * (cos(alfa)*(I(3)-I(2)/(1+P10*sin(L0)+P20*cos(L0)))...
%                 +sin(alfa)*I_test)
%             keyboard
%             
            
%             t11=3*a0^(7/2)*B^5/(muadim^(3/2))*cos(beta)*(cos(alfa)*(I(3)-I(2)/(1+P10*sin(L0)+P20*cos(L0)))-2*sin(alfa)/B*(I(8)+B/2*((I0(1)*I(2)-Ibas(1)*I(2)+corr)))) -...
%                 sqrt(a0^3/muadim)*B*(2*B^2*I(7)+3*P10*I(2))*P11*0 -...
%                 sqrt(a0^3/muadim)*B*(2*B^2*I(6)+3*P20*I(2))*P21*0;
        else
            t11=a0^(7/2)/(muadim^(3/2))*(3*cos(beta)*(sin(alfa)*((L(i)^2-L0^2)/2-(L(i)-L0)*L0)) -...
                (2*(I(4)*cos(L0)-I(1)+I(5)*sin(L0))*cos(alfa)+...
                4*(-I(4)*sin(L0)+I(5)*cos(L0))*sin(alfa))*cos(beta));
        end
%         figure(99)
%         plot(L(i),t11,'.',L(i),t11n,'o',L(i),t11-t11n,'d')

        t(i)=t00+epsilon*t11;
    end
    
%     if i == 1
%         Equin(:,1)=[a0 P10 P20 Q10 Q20 L0+(Equin0(6)-L0)]';
%     end
    % Calculate new value for the Equinoctial elements
    Equin(:,i)=[a0 P10 P20 Q10 Q20 L(i)+(Equin0(6)-L0)]'+epsilon*[a11 P11 P21 Q11 Q21 0]';
end

return

