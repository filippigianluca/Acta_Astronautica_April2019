function [Equin,t]=AnEquin_forward_m_r2(L,Equin0,epsilon,alfa,beta,muadim)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_L0 = cos(L0);
s_L0 = sin(L0);
% The following line compute e cos(theta0) where
% theta = L - (Omega + omega)
e0_c_th0 = P10*s_L0+P20*c_L0;
% The following line compute e sin(theta0) where
% theta = L - (Omega + omega)
e0_s_th0 = P20*s_L0-P10*c_L0;
% Initial true anomaly theta0
th0 = atan2(e0_s_th0, e0_c_th0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


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
dI = [0; 0; 0; Limits4_m(P10,P20)];


%% Begin cycle

Ibas = zeros(43,1);

corr = 0;

e_flag = e0 >= 1e-8;

% Limits of the analytic integrals for the discontinuity
dI = [0; 0; 0; Limits4_m(P10,P20)];

% Analytic integrals are discontinuous at L = pi.
% The following code solves the problem.
% Correction is needed only for elliptical orbits - no idea why
% Ibas will be added to I that computes the analytic integrals as
% reported in the paper.
if e_flag
    
    n_corr = 1;
    
%     if L>L0 % Case of increasing L (forward propagation)
%         while L > (n_corr*pi)
%             Ibas(4:end) = Ibas(4:end) + dI(4:end);
%             n_corr = n_corr + 2;
%             corr   = corr + dI(24)*Ic(2);
%             
%             % What is the purpose of the following line? Ic seems to be
%             % never called...
%             Ic(2)  = Ic(2) + dI(25);
%         end
%         
%     else % Case of decreasing L (backward propagation)
%         while L < -(n_corr*pi)
%             Ibas(4:end) = Ibas(4:end)-dI(4:end);
%             n_corr = n_corr+2;
%             corr = corr-dI(24)*Ic(2);
%             
%             % What is the purpose of the following line? Ic seems to be
%             % never called...
%             Ic(2) = Ic(2)-dI(25);
%         end
%     end
%     
end

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
%     I=Ibas+Integrals_m(L0,L(i),P10,P20);

    
    I =Ibas +  [0; 0; 0; Integrals4_m_r2(L0, L(i), P10, P20, 0, 0, 0, e_flag)];       
    I(24);
%    if i == 49
%        keyboard
%    end
    Phi = 1 + P10 * sin(L(i)) + P20 * cos(L(i));
    Phi0 = 1 + P10 * sin(L0) + P20 * cos(L0);
    
    
    
    if e_flag
        
        Is1 = P10 / e0^2 * (L(i) - L0) - P20 / e0^2 * log(Phi/Phi0) + ...
            -P10 / e0^2 * I(24);
        
        fun = @(x)(sin(x) ./ (1+P10*sin(x)+P20*cos(x)));
        a = L0;
        b = L(i);
%         Is1=integral(fun, a, b,'AbsTol',1e-12);
        
        Ic1 = P20 / e0^2 * (L(i) - L0) + P10 / e0^2 * log(Phi/Phi0) + ...
            -P20 / e0^2 * I(24);
        
        fun = @(x)(cos(x) ./ (1+P10*sin(x)+P20*cos(x)));
        a = L0;
        b = L(i);
%         Ic1=integral(fun, a, b,'AbsTol',1e-12);
    else
        
        Is1 = cos(L0) - cos(L(i));
        Ic1 = sin(L(i)) - sin(L0);
    end
    
%     if Ic1 < 0
%         L(i)
%         
%         keyboard
%     end
    
    a11  =  2 * a0  / (muadim * B^2) * cos(beta) * ...
        ( (P20 * cos(alfa) + P10 * sin(alfa)) * (cos(L0) - cos(L(i))) + ...
        - (P10 * cos(alfa) - P20 * sin(alfa)) * (sin(L(i)) - sin(L0)) + ...
        sin(alfa) * (L(i) - L0) );
    
    P11 = 1 / muadim * ( - cos(beta) * cos(alfa) * (sin(L(i)) - sin(L0)) + ...
        sin(alfa) * cos(beta) * (cos(L0) - cos(L(i))) + ...
        (sin(alfa) * cos(beta) + P20 * Q20 * sin(beta) ) * Is1 ...
        - Q10 * P20 * sin(beta)  * Ic1 + ...
        P10 * sin(alfa) * cos(beta) * I(24));
    
    P21 = 1 / muadim * ( cos(beta) * cos(alfa) * (cos(L0) - cos(L(i))) + ...
        sin(alfa) * cos(beta) * (sin(L(i)) - sin(L0)) + ...
        (sin(alfa) * cos(beta) + P10 * Q10 * sin(beta) ) * Ic1 ...
        - Q20 * P10 * sin(beta)  * Is1 + ...
        P20 * sin(alfa) * cos(beta) * I(24));
    
    Q11 = 1 / (2 * muadim) * (1 + Q10^2 + Q20^2) * ...
        sin(beta) * Is1;
    
    Q21 = 1 / (2 * muadim) * (1 + Q10^2 + Q20^2) * ...
        sin(beta) * Ic1;
    
    if nargout>1
        t00=h0^3*I(25);
        
        n_nodes = max([2 round((L(i)-L0)/(2*pi)*6)]);
%         n_nodes=4;
        Ls = L0 + (L(i)-L0) * linspace(0, 1, n_nodes+1);
        w  = (L(i)-L0) * [ones(n_nodes-1,1);0.5] / (0.5+n_nodes);
        dfs = dItJ5_dL(Ls(2:end), L0, a0, P10, P20, Q10, Q20, th0, ...
            a11, P11, P21);
        t11 = B * sqrt(a0 / muadim) * B * quadx_FZ(dfs,w);
        
%         [x,w]=lgwt(40,L0,L(i));
%         dfs = dItJ5_dL(x, L0, a0, P10, P20, Q10, Q20, th0, ...
%             a11, P11, P21);
%         t112 = B * sqrt(a0 / muadim) * B * quadx_FZ(dfs,w);
        
        
%         T11_TEWST = integral(@(x)function_quadrature(x, a0, ...
%         P10, P20, L0, epsilon, alfa, beta, muadim), ...
%             L0, L(i),'AbsTol',1e-12);

        t(i)=t00+epsilon*t11;
    end
    

    % Calculate new value for the Equinoctial elements
    Equin(:,i)=[a0 P10 P20 Q10 Q20 L(i)+(Equin0(6)-L0)]'+...
        epsilon*[a11 P11 P21 Q11 Q21 0]';
end

return

