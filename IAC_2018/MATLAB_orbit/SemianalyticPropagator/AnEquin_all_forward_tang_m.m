function [Equin,t,I_out] = AnEquin_all_forward_tang_m(L_values, Equin0, epsilon_t, beta_t, ...
                                                       epsilon_rth, alpha_rth, beta_rth, ...
                                                       epsilon_In, alpha_In, beta_In, ...
                                                       muadim, geopotential, third_body, ...
                                                       R, drag, Earth_flat, ill_flag, ...
                                                       flag_eps_r2)
% keyboard
% disp('int')
% Function which implements the analytical formulas for the perturbative
% expansion of Keplerian motion.

% INPUTS:
% L: row vector containing the inital value for longitude L0 and the points
%    in which the Equinoctial elements should be evaluated.
% Equin0: column vector containing the initial state (at L0) in Equinoctial
%         elements.
% epsilon_t, beta_t: modulus and elevation (in radians) of the perturbing
%                    constant tangential acceleration vector 
% epsilon_rth,alpha_rth,beta_rth: modulus, azimuth and elevation (in radians) of the
%                                 perturbing acceleration vector, in the radial, 
%                                 transverse, normal, reference system.
% epsilon_In, alpha_In, beta_In: modulus, azimuth and elevation (in radians) of the
%                                perturbing acceleration vector, costant in
%                                the inertial reference frame.
% muadim: gravitational constant of the central body, with dimensions
%         consistent with the other quantities.
% J2:
% R: Earth radius
% ill_flag: 
%
% OUTPUTS:
% Equin: Matrix containing the column vectors in Equinoctial elements for
%   each value of longitude contained in L.
% t: vector of times corresponding to the equinoctial elements.
%
% Calls "Integrals4" and "Limits4" and also others....


% REFERENCE:
% Zuiani, "Extended Analytical Formulas for the Perturbed Keplerian Motion 
% Under a Constant Control Acceleration"
%
%   Federico Zuiani 2014
%   Marilena Di Carlo, 2015-2016
%   marilena.di-carlo@strath.ac.uk

%% Initialization

if nargin < 18
    flag_eps_r2 = 0;
end

% 
% if isa( epsilon_rth_v, 'double')
%     epsilon_rth = epsilon_rth_v;
% elseif  isa( epsilon_rth_v, 'struct')
%     epsilon_rth = epsilon_rth_v.eps;
%     R_ref = epsilon_rth_v.R_ref;
% end

% Where are we going to compute the equinoctial elements? If the first
% element in L coincide with the true longitude of Equin0 and if L has more
% than one element, equinocital elements are computed considering the second
% defined true longitude in the vector L.
n = length(L_values);
nn = n - 1;

if n == 1 
    L = L_values;
    nn = n;
end

for i = 1 : nn

% if (L(1) == Equin0(6)) && (length(L)>1)
%     L = L(2);
% end
if (L_values(1) == Equin0(6)) && (length(L_values)>1)
    L = L_values(i+1);
end


% Initial true longitude [rad]
L0 = Equin0(6);

% keyboard
% If the difference between the final and initial true longitude is big
% enough....
% if L(end)-L0> 1e-8
    if abs(L(end)-L0)> 1e-8
    
    % If the number of input to the function is lower than 14, that is
    % ill_flag was not defined, automatically define ill_flag equal to 1
    if nargin < 14
       ill_flag = 1; 
    end
    
    % Initial Equinoctial Elements
    a0  = Equin0(1);
    P10 = Equin0(2);
    P20 = Equin0(3);
    Q10 = Equin0(4);
    Q20 = Equin0(5);
    % l0=Equin0(6);
    
    % Initial eccentricity
    e0 = sqrt(P10^2 + P20^2);
    
    % e_flag is 1 if e0 is greater than approximately zero
    e_flag = e0 >= 1e-8;
    
    % Parameter B (equation [3] in Reference Zuiani)
    B = sqrt(1-P10^2-P20^2);
    
    % Parameter p (initial value)
    p0 = a0*B^2;
    
    % Angular momentum (initial value)
    h0 = sqrt(muadim * p0);
    
    % Computation of the initial value of Omom = (Omega + omega)
    % Omom is the sum of RAAN and perigee argument.
    % If eccentricity is zero P10 and P20 are zero so Omom can not be
    % computed
    if e_flag
%         Omom = atan2(P10,P20);
        Omom = atan2(real(P10), real(P20));
    else
        Omom = 0;
    end
    
    % Equation [34] (J2 perturbation section)
    G = 1 + Q10^2 + Q20^2;
    
    % ?
    M = G-2;
    
    % =====================================================================
    % - Constant Acceleration in the r-t-h Frame
    % - Constant Tangential Acceleration
    % ===================================================================== 
    % Coefficient for the expression of the semimajor axis, Equations 17
    % and 30
    ka = 2*a0^3*B^2/muadim;
    
    % Coefficients for the expression of P1 and P2, Equations 17 and 30
    kP = (a0*B^2)^2/muadim;
    
    % Coefficient for the expression of Q1 and Q2, Equations 17
    kQ = kP*G/2;
    
    % Coefficient for the time in equation 19
    kt = sqrt(a0^7/muadim^3) * B^5;
  

    
    % Coefficienti per J2 probabilmente, ma non ne sono sicura
   J2 = geopotential.J2;
   J3 = geopotential.J3;
   J4 = geopotential.J4;
   J5 = geopotential.J5;
   
    k1=3*J2*R^2/(B^6*G^2*a0);
    k2=3*J2*R^2/(B^4*G^2*a0^2);
    k3=3*J2*R^2/(B^4*G*a0^2);
    
    % Coefficient for J3
    kJ3 = J3 * R^3 / (a0^3);
    
    % Coefficient for J4 
    kJ4 = J4 * R^4 / (a0^4);
    
    % Coefficient for J5
    kJ5 = J5 * R^5 /(a0^5);
    
    

    % =====================================================================
    % Warning if trajectory has eccentricity greater than 1 or if
    % periapsis is lower than Earth radius
    % =====================================================================
    if e0>1
        warning('Cannot treat open trajectories')
    end
    
    if a0*(1-e0)<=R
       % warning('Initial orbit''s pericentre is below central body''s surface')
    end

    
    

    % Non capisco. Il segente codice sembra riportare L0 tra -pi e pi,
    % invece che tra [0,2*pi) come dice di fare.....
    % Perche' calcola gli integrali usando  pi come valore per L?!?!?
    
    % L0 should be comprised between [0,2*pi)
    if L>L0
        while L0>=pi
            L0=L0-2*pi;
        end
        while L0<-pi
            L0=L0+2*pi;
        end
        
        Ic=Integrals_m(L0,pi,P10,P20);
        
    else
        while L0>pi
            L0=L0-2*pi;
        end
        while L0<=-pi
            L0=L0+2*pi;
        end
        Ic=Integrals_m(L0,-pi,P10,P20);
    end
    
    % Adjust the rest of L accordingly
    L=L-(Equin0(6)-L0);
    
    
    % =====================================================================
    % Useful precomputation
    % =====================================================================
    c_L0 = cos(L0);
    s_L0 = sin(L0);
    
    %Compute th and th0:
    % ---------------------------------------------------------------------
    % Commented by Federico:
    % th=L-Omom;
    % th0=L0-Omom;
    % ---------------------------------------------------------------------
    
    % The following line compute e cos(theta0) where 
    % theta = L - (Omega + omega)
    e0_c_th0 = P10*s_L0+P20*c_L0;
    % The following line compute e sin(theta0) where
    % theta = L - (Omega + omega)
    e0_s_th0 = P20*s_L0-P10*c_L0;
    
    % Initial true anomaly theta0
    th0 = atan2(e0_s_th0, e0_c_th0);
    
    % ---------------------------------------------------------------------
    % Commented by Federico:
    % e0_c_th=P10*sin(L)+P20*cos(L);
    % e0_s_th=P20*sin(L)-P10*cos(L);
    % ---------------------------------------------------------------------
    
    % Since (Omega + omega) remains always constant
    % L - theta = L0 - theta0
    % and therefore, the true anomaly theta is:
    th = th0 + (L-L0);
    
        
    % ---------------------------------------------------------------------
    % Commented by Federico:
    % th=atan2(e0_s_th,e0_c_th);
    % if L>L0
    %     if th<th0
    %         th=th+2*pi;
    %     end
    %     th=th+floor((L-L0)/(2*pi))*2*pi;
    % else
    %     if th>th0
    %         th=th-2*pi;
    %     end
    %     th=th-floor((L-L0)/(2*pi))*2*pi;
    % end
    % ---------------------------------------------------------------------
    
    %% Accelerations
    
    % r-t-h acceleration ????
    cbeta_t=cos(beta_t);
    sbeta_t=sin(beta_t);
    
    % r-t-h acceleration
    calpha_rth = cos(alpha_rth);
    salpha_rth = sin(alpha_rth);
    cbeta_rth  = cos(beta_rth);
    sbeta_rth  = sin(beta_rth);
    
    % Inertial acceleration orientation
    gamma0 = L0 + alpha_In;        % Equation [22]
    cgamma0 = cos(gamma0);
    sgamma0 = sin(gamma0);
    cbeta0  = cos(beta_In);
    sbeta0  = sin(beta_In);
    
    % ---------------------------------------------------------------------
    % Commented by Federico:
    % Calculation of limits for the discontinuity at L=pi+2*k*pi, k=1,2,3...
    % if epsilon_t
    %     errtol=1e-5;
    %     EllE=lellipe3(pi/2,4*e0/(1+e0)^2,errtol);%ellipticE(4*e0/(1+e0)^2);
    %     EllK=lellipf3(pi/2,4*e0/(1+e0)^2,errtol);%ellipticK(4*e0/(1+e0)^2);
    %     DIPa=(EllE-(1+e0^2)/(1+e0)^2*EllK)/(e0*(1-e0));
    %     dI=[2*(EllE/(1-e0)+EllK/(1+e0));2*P10*DIPa/e0;2*P20*DIPa/e0;Limits4_m(P10,P20)];
    % else
    % end
    % ---------------------------------------------------------------------
    


    
    %% Begin cycle
    

    
%     % Initialize equinoctial elements and time
%     Equin = NaN*ones(6,1);
%     t = NaN;
    
    % Cosa fa qui?!?!?!?!
    % Adjust the initial value for multiplier n_corr which accounts for the
    % discontinuity.
    
    % =====================================================================
    % Discontinuity correction
    % =====================================================================
    Ibas = zeros(43,1);
    corr = 0;
    
    % Limits of the analytic integrals for the discontinuity
    dI = [0; 0; 0; Limits4_m(P10,P20)];
    
    % Analytic integrals are discontinuous at L = pi. 
    % The following code solves the problem.
    % Correction is needed only for elliptical orbits - no idea why
    % Ibas will be added to I that computes the analytic integrals as
    % reported in the paper. 
    if e_flag 
        
        n_corr = 1;
        
        if L>L0 % Case of increasing L (forward propagation)
            while L > (n_corr*pi)
                Ibas(4:end) = Ibas(4:end) + dI(4:end);
                n_corr = n_corr + 2;
                corr   = corr + dI(24)*Ic(2);
                
                % What is the purpose of the following line? Ic seems to be
                % never called...
                Ic(2)  = Ic(2) + dI(25);
            end
            
        else % Case of decreasing L (backward propagation)
            while L < -(n_corr*pi)
                Ibas(4:end) = Ibas(4:end)-dI(4:end);
                n_corr = n_corr+2;
                corr = corr-dI(24)*Ic(2);
                
                % What is the purpose of the following line? Ic seems to be
                % never called...
                Ic(2) = Ic(2)-dI(25);
            end
        end
        
    end
    
    %% Analytic integrals
    
    % =====================================================================
    % WITH DRAG
    % =====================================================================
    if drag.CD
        
%         n_nodes = 1;
        n_nodes = max([2 round((L-L0)/(2*pi)*6)]);
% n_nodes = 20;
        Ls = L0 + (L-L0) * linspace(0, 1, n_nodes+1);
%         Ls = L;
        w  = (L-L0) * [ones(n_nodes-1,1);0.5] / (0.5+n_nodes);
        
        % The time has to be numerically integrated using quadrature rule.
        % To do so, we need values of the analytical itnegrals al the
        % points of the Gaussian quadrature, that's why the analytical
        % integral for the drag are computed at points Ls
        
        % OLD IMPLEMENTATION - RHO CONSTANT OVER THE ORBIT
%         [IDrag,dIDrag] = Integrals_Drag(L0, Ls(2:end), P10, P20, e_flag);

        % NEW IMPLEMENTATION - RHO DEFINED AS INTERPOLATION OF CHEBYCHEV
        % POLYNOMIALS
  
        
%         [IDrag,dIDrag] = Integrals_Drag3_2(L0, Ls(2:end), P10, P20, a0, drag, R, e_flag); 
       
        [IDrag,dIDrag] = Integrals_Drag3_2(L0, Ls(2:end), P10, P20, Q10, Q20, a0, drag, R, J2, e_flag, Earth_flat); 
        
        

        % -----------------------------------------------------------------
        % Secondo me questa parte non serve perche' i limiti destro e
        % sinistro per gli integrali del drag. Lo tengo finche' Federico
        % non conferma che e' giusto
        n_corr_th = 0;
        
        if th > th0 % Case of increasing th (forward propagation)
            while th>((1+n_corr_th*2)*pi)
                n_corr_th=n_corr_th+1;
            end
        else % Case of decreasing th (backward propagation)
            while th<((-1+n_corr_th*2)*pi)
                n_corr_th=n_corr_th-1;
            end
        end

        IDrag = n_corr_th * dIDrag + IDrag(:,end);
        % -----------------------------------------------------------------

        % Time integration for the drag acceleration case
        dfs_Drag = dIt_Drag_dL(Ls(2:end), L0, P10, P20, th0, IDrag, dIDrag(1:3));
%         dfs_Drag = dIt_Drag_dL(L, L0, P10, P20, th0, IDrag, dIDrag(1:3));

        It1_Drag = sqrt(a0  / muadim) * B * a0 * quadx_FZ(dfs_Drag,w);
    end
    
    
    
    
    
    % =====================================================================
    % WITH CONSTANT TANGENTIAL ACCELERATION
    % =====================================================================
    % Analytic integrals in the case in which there is the constant
    % tangential acceleration (elliptic integrals needs to be computed)
    if epsilon_t
        
        % Create vector with the same dimension of Ibas but with zero
        % elements
        I = Ibas*0;
        
        % -----------------------------------------------------------------
        % Commented by Federico:
        %     if th(end)>th0
        %         k=(th>pi).*(1+floor((th-pi)/(2*pi)));
        %     else
        %         k=(th<-pi).*(-1+ceil((th+pi)/(2*pi)));
        %     end
        % -----------------------------------------------------------------
        
        
        n_corr_th = 0;
        
        if th > th0 % Case of increasing th (forward propagation)
            while th>((1+n_corr_th*2)*pi)
                n_corr_th=n_corr_th+1;
            end
        else % Case of decreasing th (backward propagation)
            while th<((-1+n_corr_th*2)*pi)
                n_corr_th=n_corr_th-1;
            end
        end
        
        % -----------------------------------------------------------------
        % Commented by Federico:
        % Time integral compuation and precomputations of all the elliptic integrals needed.
        % Version A: Nodes taken from Gauss-Legendre rule
        %         if abs((L-L0)-2*pi)<1e-6
        %             x=[0.980144928248768;0.898333238706813;0.762766204958165;0.591717321247825;...
        %                 0.408282678752175;0.237233795041836;0.101666761293187;0.019855071751232];
        %             w=[0.050614268145188;0.111190517226687;0.156853322938944;0.181341891689181;...
        %                 0.181341891689181;0.156853322938944;0.111190517226687;0.050614268145188];
        %             % x=[0.990780317123360;0.952058628185237;0.884951337097152;0.793658977143309;0.683915749499090;0.562616704255734;...
        %             %     0.437383295744266;0.316084250500910;0.206341022856691;0.115048662902848;0.047941371814763;0.009219682876640];
        %             % w=[0.023587668193256;0.053469662997659;0.080039164271673;0.101583713361533;0.116746268269177;0.124573522906701;...
        %             %     0.124573522906701;0.116746268269177;0.101583713361533;0.080039164271673;0.053469662997659;0.023587668193256];
        %             Ls=L0+(L-L0)*x;
        %             w=w*(L-L0);
        %         else
        %             [Ls,w]=lgwt((n_corr_th+1)*6,L0,L);
        %         end
        %         [Is,dI(1:3)]=Integrals_tang_m(L0,[Ls.' L],P10,P20,e_flag);
        %         dfs=dIt_dL(Ls,L0,P10,P20,th0,Is(:,1:end-1),dI(1:3));
        % -----------------------------------------------------------------
        
        % Version B: Equally spaced nodes (less expensive)
        n_nodes = max([2 round((L-L0)/(2*pi)*6)]);
        Ls = L0 + (L-L0) * linspace(0, 1, n_nodes+1);
        w  = (L-L0) * [ones(n_nodes-1,1);0.5] / (0.5+n_nodes);
        [Is, dI(1:3)] = Integrals_tang_m(L0, Ls(2:end), P10, P20, e_flag);
        
        % The output of Is is, for each value of Ls, is:
        % Is = [Ia Ip1 Ip2] where Ip1 and Ip2 are not the one in the paper
        % but:
        % IP1 = s_Omom*IPa + c_Omom*IPb;
        % IP2 = c_Omom*IPa-s_Omom*IPb;
        % and IPa and IPb are what in the paper are called Ip1 and Ip2
        
        % dfs is the function to be integrated to obtain t1 as in Equation
        % 31. The value of the integrand f for each value of Ls is stored
        % in dfs. These values are later used for the integration with
        % quadrature rule using the weight w.
        dfs = dIt_dL(Ls(2:end), L0, P10, P20, th0, Is, dI(1:3));
        It = quadx_FZ(dfs,w);
  
        % Vector of analytic integrals: the first three elements are the
        % analytic integrals for the case of constant tangential thrust.
        % The remaining are the one computed using Integrals4_m and
        % corrected for the discontinuity in pi using Ibas

        I(1:3)   = n_corr_th * dI(1:3) + Is(:,end);      
        I(4:end) = Ibas(4:end) + Integrals4_m(L0, L, P10, P20, epsilon_rth~=0, epsilon_In~=0, J2~=0, e_flag);
        
    % =====================================================================
    % WITHOUT CONSTANT TANGENTIAL ACCELERATION
    % =====================================================================
    % Analytic integrals in the case in which there is NOT the constant
    % tangential acceleration
    else
    
        I = Ibas + [0; 0; 0; Integrals4_m(L0, L, P10, P20, epsilon_rth~=0, epsilon_In~=0, J2~=0, e_flag)];       
        I_out(:,i) = I;
%         The elements in I are:
%         I = [0; 0; 0; [IJc1;IJs1];...
%             [IJc2;IJc1s1;IJs2];...
%             [IJc3;IJc2s1;IJc1s2;IJs3];...
%             [IJc4;IJc3s1;IJc2s2;IJc1s3;IJs4];...
%             [IJc5;IJc4s1;IJc3s2;IJc2s3;IJc1s4;IJs5];...
%             [I11;I12;I13;Ic2;Is2;Ic3;Is3;Ia2];     from 24 to 31
%             [Icc3;Ics3;Iss3;I3c3;I2c1s3;I1c2s3;I3s3;...   from 32
%              I4c3;I3c1s3;I2c2s3;I1c3s3;I4s3]];
%              

    end
    
    
    %% Constant Tangential Acceleration
    % Equations [30] - ma non tornano, non capisco. Come fa ad esserci
    % beta?!?!?!?!?
    % Devono essere quelle tangenziali, perche` altrimenti dove sono?!?!?
    
    % Calculate first order variation of Equinoctial elements for constant
    % tangential thrust
    % Equation [30] have been generated using beta=0, so how come that now
    % we have beta?!?!?!
    % I(1) = Ia;
    % I(2) is the integral in equation [28b], that is:
    % sin(Omega+omega) Ip1 + cos(Omega+omega) Ip2
    % I(3) is the integral in equation [28b], that is:
    % cos(Omega+omega) Ip1 -sin (Omega+omega) Ip2
    a1t  = ka * cbeta_t * I(1);
    P11t = kP * 2 * cbeta_t * I(2);
    P21t = kP * 2 * cbeta_t * I(3);
    Q11t = kQ * sbeta_t * I(30);
    Q21t = kQ * sbeta_t * I(29);
    
    % Secondo me le equazioni dovrebbero essere:
    %     a1t  = ka * I(1);
    %     P11t = kP * 2 *  I(2);
    %     P21t = kP * 2 *  I(3);
    %     Q11t = Q10;
    %     Q21t = Q20;
    
    %% r-t-h thrust - verificati 
    
    % Considering equations [17], we will need:
    % I11 = I(24);
    % Ic2 = I(27);
    % I13 = I(26);
    % Is3 = I(30);
    % Is2 = I(28);
    % Ic3 = I(29);
    
    if epsilon_rth && ~flag_eps_r2
        a1rth  = ka * cbeta_rth * (calpha_rth*(P20*I(28)-P10*I(27)) + salpha_rth * I(24));   
        P11rth = kP *( cbeta_rth * (-calpha_rth*I(27) + salpha_rth*(P10*I(26)+I(30)+I(28)) ) + sbeta_rth * P20*(-Q10*I(29)+Q20*I(30)));
        P21rth = kP *( cbeta_rth * ( calpha_rth*I(28) + salpha_rth*(P20*I(26)+I(29)+I(27)) ) + sbeta_rth * P10*(Q10*I(29)-Q20*I(30)));
        Q11rth = kQ * sbeta_rth * I(30);
        Q21rth = kQ * sbeta_rth * I(29);
    else
        a1rth  = 0;
        P11rth = 0;
        P21rth = 0;
        Q11rth = 0;
        Q21rth = 0;
    end
    
    
     %% r-t-h thrust - epsilon as 1/r2
    
    % Considering equations [17], we will need:
    % I11 = I(24);
    % Ic2 = I(27);
    % I13 = I(26);
    % Is3 = I(30);
    % Is2 = I(28);
    % Ic3 = I(29);
    
    if epsilon_rth && flag_eps_r2
        
        Phi = 1 + P10 * sin(L) + P20 * cos(L);
        Phi0 = 1 + P10 * sin(L0) + P20 * cos(L0);
        
        if e_flag
            
            Is1 = P10 / e0^2 * (L - L0) - P20 / e0^2 * log(Phi/Phi0) + ...
                -P10 / e0^2 * I(24);
            
            Ic1 = P20 / e0^2 * (L - L0) + P10 / e0^2 * log(Phi/Phi0) + ...
                -P20 / e0^2 * I(24);
        else
            
            Is1 = cos(L0) - cos(L);
            Ic1 = sin(L) - sin(L0);
        end
        
        a1rth_r2  =  2 * a0  / (muadim * B^2) * cbeta_rth * ...
            ( (P20 * calpha_rth + P10 * salpha_rth) * (cos(L0) - cos(L)) + ...
            - (P10 * calpha_rth - P20 * salpha_rth) * (sin(L) - sin(L0)) + ...
            salpha_rth * (L - L0) );
        
        P11rth_r2 = 1 / muadim * ( - cbeta_rth * calpha_rth * (sin(L) - sin(L0)) + ...
            salpha_rth * cbeta_rth * (cos(L0) - cos(L)) + ...
            (salpha_rth * cbeta_rth + P20 * Q20 * sbeta_rth ) * Is1 ...
            - Q10 * P20 * sbeta_rth  * Ic1 + ...
            P10 * salpha_rth * cbeta_rth * I(24));
        
        P21rth_r2 = 1 / muadim * ( cbeta_rth * calpha_rth * (cos(L0) - cos(L)) + ...
            salpha_rth * cbeta_rth * (sin(L) - sin(L0)) + ...
            (salpha_rth * cbeta_rth + P10 * Q10 * sbeta_rth ) * Ic1 ...
            - Q20 * P10 * sbeta_rth  * Is1 + ...
            P20 * salpha_rth * cbeta_rth * I(24));
        
        Q11rth_r2 = 1 / (2 * muadim) * (1 + Q10^2 + Q20^2) * ...
            sbeta_rth * Is1;
        
        Q21rth_r2 = 1 / (2 * muadim) * (1 + Q10^2 + Q20^2) * ...
            sbeta_rth * Ic1;
    else
        a1rth_r2  = 0;
        P11rth_r2 = 0;
        P21rth_r2 = 0;
        Q11rth_r2 = 0;
        Q21rth_r2 = 0;
    end
    
    %% Inertial thrust - verificati
    
    % Considering equations [23], we will need:
    % I12 = I(25);
    % Is2 = I(28);
    % Ic2 = I(27);
    % Is3 = I(30);
    % Ic3 = I(29);
    % I2s3 (Iss3) = I(34);
    % I1c1s3 (Ics3) = I(33);
    % I2c3 (Icc3) = I(32);
    
    if epsilon_In
%         keyboard
        a1In  = ka * cbeta0 * (sgamma0*(P20*I(25) + I(27)) - cgamma0*(P10*I(25) + I(28)));
        P11In = kP * (cbeta0 * ( (P10*I(29) + I(33))*sgamma0 - (I(25) + I(34) + P10*I(30))*cgamma0) - sbeta0*(P20*Q10*I(29) - P20*Q20*I(30)));
        P21In = kP * (cbeta0 * ( (I(25) + I(32) + P20*I(29))*sgamma0 - (I(33) + P20*I(30))*cgamma0) + sbeta0*(P10*Q10*I(29) - P10*Q20*I(30)));
        Q11In = kQ * sbeta0 * I(30);
        Q21In = kQ * sbeta0 * I(29);
    else
        a1In  = 0;
        P11In = 0;
        P21In = 0;
        Q11In = 0;
        Q21In = 0;
    end
    
    %% J2 perturbation 
    if J2
        
        % Terms in the equation for the semimajor axis:
        % I(18) -> IJc5 = 5/8 * (sL-sL0) + 5/48*(s3L-s3L0) + 1/80*(s5L-s5L0);
        % I(19) -> IJc4s1 = -(cL^5 - cL0^5)/5;
        % I(13) -> IJc4   = 3/8*(L-L0) + 0.25*(s2L-s2L0) + 1/32*(s4L-s4L0);
        % I(20) -> IJc3s2 = 1/8*(sL-sL0) - 1/48*(s3L-s3L0) - 1/80*(s5L-s5L0);
        % I(14) -> IJc3s1 = -(cL^4 - cL0^4)/4;
        % I(9)  -> IJc3   = 0.75*(sL-sL0) + 1/12*(s3L-s3L0);
        % I(21) -> IJc2s3 = -1/8*(cL-cL0) - 1/48*(c3L-c3L0) + 1/80*(c5L-c5L0);
        % I(15) -> IJc2s2 = (L - L0)/8 - (s4L - s4L0)/32;
        % I(10) -> IJc2s1 = -(cL^3 - cL0^3)/3;
        % I(6)  -> IJc2   = 0.5*(L - L0) + 0.25*(s2L-s2L0);
        % I(22) -> IJc1s4 = (sL^5 - sL0^5)/5;
        % I(16) -> IJc1s3 = (sL^4 - sL0^4)/4;
        % I(11) -> IJc1s2 = (sL^3 - sL0^3)/3;
        % I(7)  -> IJc1s1 = -0.5*(cL^2-cL0^2);
        % I(4)  -> IJc1 = sL - sL0;
        % I(23) -> IJs5   = -5/8*(cL-cL0) + 5/48*(c3L-c3L0) - 1/80*(c5L-c5L0);
        % I(17) -> IJs4   = 3/8*(L-L0)-0.25*(s2L-s2L0) + 1/32*(s4L-s4L0);
        % I(12) -> IJs3   = -0.75*(cL-cL0) + 1/12*(c3L-c3L0);
        % I(8)  -> IJs2   = 0.5*(L - L0) - 0.25*(s2L - s2L0);
        % I(5)  -> IJs1 = cL0 - cL;
        a1J = k1 * ( (8*Q20*P20^3*Q10 - 12*P10*P20^2*Q10^2) * I(18) +...                                                      
            (- 24*P10^2*P20*Q10^2 + 48*P10*P20^2*Q10*Q20 + 20*P20^3*Q10^2 - 8*P20^3*Q20^2) * I(19) +...                       
            (- 24*P10*P20*Q10^2 + 24*P20^2*Q10*Q20) * I(13) +...
            (- 12*P10^3*Q10^2 + 72*P10^2*P20*Q10*Q20 + 48*P10*P20^2*Q10^2 - 36*P10*P20^2*Q20^2 - 32*P20^3*Q10*Q20) * I(20) +...
            (- 24*P10^2*Q10^2 + 96*P10*P20*Q10*Q20 + 48*P20^2*Q10^2 - 24*P20^2*Q20^2) * I(14) +...
            (P10*G^2*P20^2 + 24*Q20*P20*Q10 - 12*P10*Q10^2) * I(9) +...
            (32*P10^3*Q10*Q20 + 36*P10^2*P20*Q10^2 - 48*P10^2*P20*Q20^2 - 72*P10*P20^2*Q10*Q20 + 12*P20^3*Q20^2) * I(21) +...
            (72*P10^2*Q10*Q20 + 72*P10*P20*Q10^2 - 72*P10*P20*Q20^2 - 72*P20^2*Q10*Q20) * I(15) +...
            (2*G^2*P10^2*P20 - G^2*P20^3 + 48*P10*Q10*Q20 + 36*P20*Q10^2 - 24*P20*Q20^2) * I(10) +...
            (2*P10*P20*G^2 + 8*Q10*Q20) * I(6) +...
            (8*P10^3*Q10^2 - 20*P10^3*Q20^2 - 48*P10^2*P20*Q10*Q20 + 24*P10*P20^2*Q20^2) * I(22) +...
            (24*P10^2*Q10^2 - 48*P10^2*Q20^2 - 96*P10*P20*Q10*Q20 + 24*P20^2*Q20^2) * I(16) +...
            (G^2*P10^3 - 2*G^2*P10*P20^2 + 24*P10*Q10^2 - 36*P10*Q20^2 - 48*P20*Q10*Q20) * I(11) +...
            (2*G^2*P10^2 - 2*G^2*P20^2 + 8*Q10^2 - 8*Q20^2) * I(7) +...
            G^2*P10 * I(4) +...
            (- 8*Q10*P10^3*Q20 + 12*P20*P10^2*Q20^2) * I(23) +...
            (- 24*Q10*P10^2*Q20 + 24*P20*P10*Q20^2) * I(17) +...
            (- P20*G^2*P10^2 - 24*Q10*P10*Q20 + 12*P20*Q20^2) * I(12) +...
            (- 2*P10*P20*G^2 - 8*Q10*Q20) * I(8) +...
            (-G^2*P20) * I(5));
        %     a1J=-1;
        
        % Terms in the equation for P1:
        % I(18) -> IJc5 = 5/8 * (sL-sL0) + 5/48*(s3L-s3L0) + 1/80*(s5L-s5L0);
        % I(19) -> IJc4s1 = -(cL^5 - cL0^5)/5;
        % I(13) -> IJc4   = 3/8*(L-L0) + 0.25*(s2L-s2L0) + 1/32*(s4L-s4L0);
        % I(20) -> IJc3s2 = 1/8*(sL-sL0) - 1/48*(s3L-s3L0) - 1/80*(s5L-s5L0);
        % I(14) -> IJc3s1 = -(cL^4 - cL0^4)/4;
        % I(9)  -> IJc3   = 0.75*(sL-sL0) + 1/12*(s3L-s3L0);
        % I(21) -> IJc2s3 = -1/8*(cL-cL0) - 1/48*(c3L-c3L0) + 1/80*(c5L-c5L0);
        % I(15) -> IJc2s2 = (L - L0)/8 - (s4L - s4L0)/32;
        % I(10) -> IJc2s1 = -(cL^3 - cL0^3)/3;
        % I(6)  -> IJc2   = 0.5*(L - L0) + 0.25*(s2L-s2L0);
        % I(22) -> IJc1s4 = (sL^5 - sL0^5)/5;
        % I(16) -> IJc1s3 = (sL^4 - sL0^4)/4;
        % I(11) -> IJc1s2 = (sL^3 - sL0^3)/3;
        % I(7)  -> IJc1s1 = -0.5*(cL^2-cL0^2);
        % I(4)  -> IJc1 = sL - sL0;
        % I(23) -> IJs5   = -5/8*(cL-cL0) + 5/48*(c3L-c3L0) - 1/80*(c5L-c5L0);
        % I(17) -> IJs4   = 3/8*(L-L0)-0.25*(s2L-s2L0) + 1/32*(s4L-s4L0);
        % I(12) -> IJs3   = -0.75*(cL-cL0) + 1/12*(c3L-c3L0);
        % I(8)  -> IJs2   = 0.5*(L - L0) - 0.25*(s2L - s2L0);
        % I(5)  -> IJs1 = cL0 - cL;
        P11J = k2*( (-6*P20^2*Q10^2) * I(18) +...
            (16*P20^2*Q10*Q20 - 12*P10*P20*Q10^2) * I(19) -...
            12*P20*Q10^2 * I(13) +...
            (- 6*P10^2*Q10^2 + 32*P10*P20*Q10*Q20 + 4*P20^2*Q10^2 - 10*P20^2*Q20^2) * I(20) +...
            (- 12*P10*Q10^2 + 36*P20*Q20*Q10) * I(14) +...
            ((G^2*P20^2)/2 + 2*P20^2*Q10^4 + 2*P20^2*Q10^2*Q20^2 - 2*P20^2*Q10^2 + 4*P10*P20*Q10*Q20 - 6*Q10^2) * I(9) +...
            (16*P10^2*Q10*Q20 + 8*P10*P20*Q10^2 - 20*P10*P20*Q20^2 - 4*P20^2*Q10*Q20) * I(21) +...
            (12*P20*Q10^2 + 36*P10*Q10*Q20 - 24*P20*Q20^2) * I(15) +...
            (G^2*P10*P20 + 4*P10^2*Q10*Q20 + 2*P10*P20*Q10^4 + 2*P10*P20*Q10^2*Q20^2 + 2*P10*P20*Q10^2 - 4*P10*P20*Q20^2 - 4*P20^2*Q10^3*Q20 - 4*P20^2*Q10*Q20^3 + 4*P20^2*Q10*Q20 + 20*Q10*Q20) * I(10) +...
            (P20*G^2 + 2*P20*Q10^4 + 2*P20*Q10^2*Q20^2 - 2*P20*Q10^2 + 4*P10*Q10*Q20) * I(6) +...
            (4*P10^2*Q10^2 - 10*P10^2*Q20^2 - 8*P20*P10*Q10*Q20) * I(22) +...
            (12*P10*Q10^2 - 12*P20*Q10*Q20 - 24*P10*Q20^2) * I(16) +...
            ((G^2*P10^2)/2 + 4*P10^2*Q10^2 - 4*P10^2*Q20^2 - 4*P10*P20*Q10^3*Q20 - 4*P10*P20*Q10*Q20^3 + 2*P20^2*Q10^2*Q20^2 + 2*P20^2*Q20^4 - 2*P20^2*Q20^2 + 8*Q10^2 - 14*Q20^2) * I(11) +...
            (P10*G^2 - 4*P20*Q10^3*Q20 + 4*P10*Q10^2 - 4*P20*Q10*Q20^3 + 4*P20*Q10*Q20 - 4*P10*Q20^2) * I(7) +...
             G^2/2 * I(4) +...
            (-4*P10^2*Q10*Q20) * I(23) +...
            (-12*P10*Q10*Q20) * I(17) +...
            (- 4*P10^2*Q10*Q20 + 2*P20*P10*Q10^2*Q20^2 + 2*P20*P10*Q20^4 - 2*P20*P10*Q20^2 - 8*Q10*Q20) * I(12) +...
            (2*P20*Q10^2*Q20^2 - 4*P10*Q10*Q20 + 2*P20*Q20^4 - 2*P20*Q20^2) * I(8));
        
        P21J = k2*((4*P20^2*Q10*Q20) * I(18) +...
            (10*P20^2*Q10^2 - 4*P20^2*Q20^2 + 8*P10*P20*Q10*Q20) * I(19) +...
            12*P20*Q10*Q20 * I(13) +...
            (4*P10^2*Q10*Q20 + 20*P10*P20*Q10^2 - 8*P10*P20*Q20^2 - 16*P20^2*Q10*Q20) * I(20) +...
            (24*P20*Q10^2 + 12*P10*Q10*Q20 - 12*P20*Q20^2) * I(14) +...
            (4*P20^2*Q10*Q20 - 2*P10*P20*Q10^4 - 2*P10*P20*Q10^2*Q20^2 + 2*P10*P20*Q10^2 + 8*Q10*Q20) * I(9) +...
            (10*P10^2*Q10^2 - 4*P10^2*Q20^2 - 32*P10*P20*Q10*Q20 + 6*P20^2*Q20^2) * I(21) +...
            (24*P10*Q10^2 - 36*P20*Q10*Q20 - 12*P10*Q20^2) * I(15) +...
            (- 2*P10^2*Q10^4 - 2*P10^2*Q10^2*Q20^2 + 2*P10^2*Q10^2 + 4*P10*P20*Q10^3*Q20 + 4*P10*P20*Q10*Q20^3 + 4*P20^2*Q10^2 - 4*P20^2*Q20^2 - (G^2*P20^2)/2 + 14*Q10^2 - 8*Q20^2) * I(10) +...
            (- 2*P10*Q10^4 - 2*P10*Q10^2*Q20^2 + 2*P10*Q10^2 + 4*P20*Q10*Q20) * I(6) +...
            (- 16*Q10*P10^2*Q20 + 12*P20*P10*Q20^2) * I(22) +...
            (12*P20*Q20^2 - 36*P10*Q10*Q20) * I(16) +...
            (- G^2*P10*P20 + 4*P10^2*Q10^3*Q20 + 4*P10^2*Q10*Q20^3 - 4*P10^2*Q10*Q20 - 2*P10*P20*Q10^2*Q20^2 + 4*P10*P20*Q10^2 - 2*P10*P20*Q20^4 - 2*P10*P20*Q20^2 - 4*P20^2*Q10*Q20 - 20*Q10*Q20) * I(11) +...
            (- P20*G^2 + 4*P10*Q10^3*Q20 + 4*P20*Q10^2 + 4*P10*Q10*Q20^3 - 4*P10*Q10*Q20 - 4*P20*Q20^2) * I(7) +...
            (6*P10^2*Q20^2) * I(23) +...
            (12*P10*Q20^2) * I(17) +...
            (- 2*P10^2*Q10^2*Q20^2 - 2*P10^2*Q20^4 + 2*P10^2*Q20^2 - (G^2*P10^2)/2 - 4*P20*P10*Q10*Q20 + 6*Q20^2) * I(12) +...
            (- P10*G^2 - 2*P10*Q10^2*Q20^2 - 4*P20*Q10*Q20 - 2*P10*Q20^4 + 2*P10*Q20^2) * I(8) +...
            (-G^2/2) * I(5));
        
        Q11J = k3*((- P20*Q10^3 - P20*Q10*Q20^2 + P20*Q10) * I(10) +...
            (- P10*Q10^3 + P20*Q10^2*Q20 - P10*Q10*Q20^2 + P10*Q10 + P20*Q20^3 - P20*Q20) * I(11) +...
            (- Q10^3 - Q10*Q20^2 + Q10) * I(7) +...
            (P10*Q10^2*Q20 + P10*Q20^3 - P10*Q20) * I(12) +...
            (Q10^2*Q20 + Q20^3 - Q20) * I(8));
        
        % Terms in the equation for Q2:
        % I(9)  -> IJc3   = 0.75*(sL-sL0) + 1/12*(s3L-s3L0);
        % I(10) -> IJc2s1 = -(cL^3 - cL0^3)/3;
        % I(6)  -> IJc2   = 0.5*(L - L0) + 0.25*(s2L-s2L0);
        % I(11) -> IJc1s2 = (sL^3 - sL0^3)/3;
        % I(7)  -> IJc1s1 = -0.5*(cL^2-cL0^2);
        Q21J = k3*((- P20*Q10^3 - P20*Q10*Q20^2 + P20*Q10) * I(9) +...
            (- P10*Q10^3 + P20*Q20^3 + P10*Q10 - P20*Q20 - P10*Q10*Q20^2 + P20*Q10^2*Q20) * I(10) +...
            (Q10 - Q10*Q20^2 - Q10^3 ) * I(6) +...
            (P10*Q10^2*Q20 + P10*Q20^3 - P10*Q20) * I(11) +...
            (Q10^2*Q20 + Q20^3 - Q20) * I(7));

    else
        a1J  = 0;
        P11J = 0;
        P21J = 0;
        Q11J = 0;
        Q21J = 0;
    end
    
    
    
    %% J3 perturbation
    
    if J3
        % Analytic integrals:
        I_J3_a1 = Integral_J3_a1(Equin0, L0, L);
        I_J3_a2 = Integral_J3_a2(Equin0, L0, L);
        I_J3_P11 = Integral_J3_P1_1(Equin0, L0, L);
        I_J3_P12 = Integral_J3_P1_2(Equin0, L0, L);
        I_J3_P13 = Integral_J3_P1_3(Equin0, L0, L);
        I_J3_P21 = Integral_J3_P2_1(Equin0, L0, L);
        I_J3_Q1  = Integral_J3_Q1(Equin0, L0, L);
        I_J3_Q2  = Integral_J3_Q2(Equin0, L0, L);
        
        % Equations for the orbital elements
        a1J3 = kJ3 * a0 / (B^8 * G) * ( 8 *  I_J3_a1 + 6 * I_J3_a2);
        P11J3 = kJ3 / (B^6 * G) * (-4 * I_J3_P11 + 3 * I_J3_P12) + ...
               - 3 * kJ3 / (2 * B^6 * G) * P20 * (1 - Q10^2 - Q20^2) * I_J3_P13;
        P21J3 = kJ3 / (B^6 * G) * I_J3_P21 + ...
               + 3 * kJ3 / (2 * B^6 * G) * P10 * (1 - Q10^2 - Q20^2) * I_J3_P13;
        Q11J3 = 3/4 * kJ3 / (B^6 * G) * (1 + Q10^2 + Q20^2) * (1 - Q10^2 - Q20^2) * I_J3_Q1;
        Q21J3 = 3/4 * kJ3 / (B^6 * G) * (1 + Q10^2 + Q20^2) * (1 - Q10^2 - Q20^2) * I_J3_Q2;
    else
        a1J3 = 0;
        P11J3 = 0;
        P21J3 = 0;
        Q11J3 = 0;
        Q21J3 = 0;
    end
    
    
    %% J4 perturbation
    
    if J4
        % Analytic integrals:
        I_J4_a1 = Integral_J4_a1(Equin0, L0, L);
        I_J4_a2 = Integral_J4_a2(Equin0, L0, L);
        I_J4_P11 = Integral_J4_P1_1(Equin0, L0, L);
        I_J4_P12 = Integral_J4_P1_2(Equin0, L0, L);
        I_J4_P21 = Integral_J4_P2_1(Equin0, L0, L);
        I_J4_Q1  = Integral_J4_Q1(Equin0, L0, L);
        I_J4_Q2  = Integral_J4_Q2(Equin0, L0, L);
        
        % Equations for the orbital elements
        a1J4 = kJ4 * a0 / (B^10) * ( 0.25 *  I_J4_a1 - 2 * I_J4_a2/G^2);
        P11J4 = kJ4 / (B^8) * (-I_J4_P11 ) + ...
               -  kJ4 / ( B^8 * G^2) * P20 * (1 - Q10^2 - Q20^2) * I_J4_P12;
        P21J4 = kJ4 / (B^8) * (I_J4_P21 ) + ...
               +  kJ4 / ( B^8 * G^2) * P10 * (1 - Q10^2 - Q20^2) * I_J4_P12;
        Q11J4 = -0.5 * kJ4 / (B^8 * G^2) * (1 + Q10^2 + Q20^2) * (1 - Q10^2 - Q20^2) * I_J4_Q1;
        Q21J4 =  -0.5 * kJ4 / (B^8 * G^2) * (1 + Q10^2 + Q20^2) * (1 - Q10^2 - Q20^2) * I_J4_Q2;
    else
        a1J4 = 0;
        P11J4 = 0;
        P21J4 = 0;
        Q11J4 = 0;
        Q21J4 = 0;
    end
    
    
    
        %% J5 perturbation
    
    if J5
        % Analytic integrals:
        I_J5_a = Integral_J5_a(Equin0, L0, L);
        I_J5_P11 = Integral_J5_P1_1(Equin0, L0, L);
        I_J5_P12 = Integral_J5_P1_2(Equin0, L0, L);
        I_J5_P21 = Integral_J5_P2_1(Equin0, L0, L);
        I_J5_Q1  = Integral_J5_Q1(Equin0, L0, L);
        I_J5_Q2  = Integral_J5_Q2(Equin0, L0, L);
        
        % Equations for the orbital elements
        a1J5 = kJ5 * a0 / (B^12 * G) * I_J5_a;
        P11J5 = - kJ5 / (2 * G * B^10) * I_J5_P11  + ...
               + kJ5 / ( 8 * B^10 * G) * P20 * (1 - Q10^2 - Q20^2) * I_J5_P12;
        P21J5 = - kJ5 / (2 * G * B^10) * (I_J5_P21 ) + ...
               -  kJ5 / ( 8 * B^10 * G) * P10 * (1 - Q10^2 - Q20^2) * I_J5_P12;
        Q11J5 = -1/16 * kJ5 / (B^10 * G) * (1 + Q10^2 + Q20^2) * (1 - Q10^2 - Q20^2) * I_J5_Q1;
        Q21J5 =  -1/16 * kJ5 / (B^10 * G) * (1 + Q10^2 + Q20^2) * (1 - Q10^2 - Q20^2) * I_J5_Q2;
    else
        a1J5 = 0;
        P11J5 = 0;
        P21J5 = 0;
        Q11J5 = 0;
        Q21J5 = 0;
    end
    
    
    %% Drag acceleration
    
    % The analytical integrals IDrag were computed at many points between
    % L0 and L because of the gaussian quadrature. Now we need the one
    % corresponding to the last value L, so IDrag(:,end) is considered in
    % the following equations
    
    if drag.CD
        
        drag.k = -  drag.CD * drag.A_m;
        drag.eps =  a0 * drag.k;
        
        a1Drag = a0 * (e0^2 * IDrag(1,end) + IDrag(2,end));
        P11Drag = 0.5 * B^2 * (sin(Omom) * (e0 * IDrag(1,end) + e0 * IDrag(4,end) + 2 * IDrag(6,end) + e0 * IDrag(7,end)) + 2 * cos(Omom) * IDrag(5,end));
        P21Drag = 0.5 * B^2 * (cos(Omom) * (e0 * IDrag(1,end) + e0 * IDrag(4,end) + 2 * IDrag(6,end) + e0 * IDrag(7,end)) - 2 * sin(Omom) * IDrag(5,end));
        Q11Drag = 0;
        Q21Drag = 0;
%         keyboard
    else
        
        drag.eps = 0;
        a1Drag = 0;
        P11Drag = 0;
        P21Drag = 0;
        Q11Drag = 0;
        Q21Drag = 0;
    end
    
    %% Third body Sun
    
    
    if third_body.flag_Sun
        
        R_3rd_Sun_vector = third_body.R_3rd_Sun_vector;
        R_3rd_Sun        = norm(R_3rd_Sun_vector);
        mu_3rd_Sun       = third_body.mu_3rd_Sun;
        
        f = 1 / G * [1 - Q10^2 + Q20^2; 2 * Q10 * Q20; -2 * Q10];
        g = 1 / G * [2 * Q10 * Q20; 1 + Q10^2 - Q20^2; 2 * Q20];
        w = 1 / G * [2 * Q10; - 2 * Q20; 1 - Q10^2 - Q20^2];
        
        alpha_Sun = dot(f, R_3rd_Sun_vector) / norm(R_3rd_Sun_vector);
        beta_Sun  = dot(g, R_3rd_Sun_vector) / norm(R_3rd_Sun_vector);
        gamma_Sun = dot(w, R_3rd_Sun_vector) / norm(R_3rd_Sun_vector);
        
        eps_3rd_Sun = (a0 / R_3rd_Sun)^3;
     
        I_3rd = I(25);
        
        [coeffP11, coeffP21, coeffQ11, coeffQ21] = coefficient_integral_3rd_body(Equin0, alpha_Sun, beta_Sun, gamma_Sun, R_3rd_Sun);
        
        a1_3rd_Sun = 0;
        P11_3rd_Sun =    B^2 / (256 * R_3rd_Sun^4) * (mu_3rd_Sun / muadim) * coeffP11 * I_3rd;
        P21_3rd_Sun = -  B^2 / (256 * R_3rd_Sun^4) * (mu_3rd_Sun / muadim) * coeffP21 * I_3rd;
        Q11_3rd_Sun =    B^2 * gamma_Sun * G / (2 * 256 * R_3rd_Sun^4) * (mu_3rd_Sun / muadim) * coeffQ11 * I_3rd;
        Q21_3rd_Sun =    B^2 * gamma_Sun * G / (2 * 256 * R_3rd_Sun^4) * (mu_3rd_Sun / muadim) * coeffQ21 * I_3rd;
        
%         keyboard
        
    else
        
        eps_3rd_Sun = 0;
        
        a1_3rd_Sun  = 0;
        P11_3rd_Sun = 0;
        P21_3rd_Sun = 0;
        Q11_3rd_Sun = 0;
        Q21_3rd_Sun = 0;
        
    end
    
    
  %% Third body Moon
    
    
    if third_body.flag_Moon
        
        R_3rd_Moon_vector = third_body.R_3rd_Moon_vector;
        R_3rd_Moon        = norm(R_3rd_Moon_vector);
        mu_3rd_Moon       = third_body.mu_3rd_Moon;
        
        f = 1 / G * [1 - Q10^2 + Q20^2; 2 * Q10 * Q20; -2 * Q10];
        g = 1 / G * [2 * Q10 * Q20; 1 + Q10^2 - Q20^2; 2 * Q20];
        w = 1 / G * [2 * Q10; - 2 * Q20; 1 - Q10^2 - Q20^2];
        
        alpha_Moon = dot(f, R_3rd_Moon_vector) / norm(R_3rd_Moon_vector);
        beta_Moon  = dot(g, R_3rd_Moon_vector) / norm(R_3rd_Moon_vector);
        gamma_Moon = dot(w, R_3rd_Moon_vector) / norm(R_3rd_Moon_vector);
        
        eps_3rd_Moon = (a0 / R_3rd_Moon)^3;
     
        I_3rd = I(25);
        
        [coeffP11, coeffP21, coeffQ11, coeffQ21] = coefficient_integral_3rd_body(Equin0, alpha_Moon, beta_Moon, gamma_Moon, R_3rd_Moon);
        
        a1_3rd_Moon = 0;
        P11_3rd_Moon =    B^2 / (256 * R_3rd_Moon^4) * (mu_3rd_Moon / muadim) * coeffP11 * I_3rd;
        P21_3rd_Moon = -  B^2 / (256 * R_3rd_Moon^4) * (mu_3rd_Moon / muadim) * coeffP21 * I_3rd;
        Q11_3rd_Moon =    B^2 * gamma_Moon * G / (2 * 256 * R_3rd_Moon^4) * (mu_3rd_Moon / muadim) * coeffQ11 * I_3rd;
        Q21_3rd_Moon =    B^2 * gamma_Moon * G / (2 * 256 * R_3rd_Moon^4) * (mu_3rd_Moon / muadim) * coeffQ21 * I_3rd;
        
%         keyboard
        
    else
        
        eps_3rd_Moon = 0;
        
        a1_3rd_Moon  = 0;
        P11_3rd_Moon = 0;
        P21_3rd_Moon = 0;
        Q11_3rd_Moon = 0;
        Q21_3rd_Moon = 0;
        
    end
% %     
    %% Very illegal thing

    % Cosa fa?
    if ill_flag
        if J2~=0
            e1=sqrt((P10+P11J)^2+(P20+P21J)^2);
            if e1~=0
                P11J=e0*(P10+P11J)/e1-(P10);
                P21J=e0*(P20+P21J)/e1-(P20);

            end
            tg_i_0=sqrt(Q10^2+Q20^2);
            tg_i_1=sqrt((Q10+Q11J)^2+(Q20+Q21J)^2);
            if tg_i_1~=0
                Q11J=tg_i_0*(Q10+Q11J)/tg_i_1-Q10;
                Q21J=tg_i_0*(Q20+Q21J)/tg_i_1-Q20;
            end
        end
    end

    
    %% Time 
    
    % If output are more than one, therefore time has to be an output of
    % the function:
    if nargout > 1
        
        t00 = h0^3 * I(25);
        
        % Modified analytic
        
        % If initial eccentricity is greater than 0
        if e_flag
            
            % WITH CONSTANT TANGENTIAL ACCELERATION
            if epsilon_t
                % ---------------------------------------------------------
                % Commented by Federico:
                %             It=integral(@(LL) dIt_dL(LL,L0,P10,P20,dI(1:3)),L0,L,'RelTol',1e-2);
                %             It=quadgk(@(LL) dIt_dL(LL,L0,P10,P20,dI(1:3)),L0,L,'MaxIntervalCount',((n_corr-1)/2+1)*8);
                %             It=quadgl_FZ(@(LL) dIt_dL(LL,L0,P10,P20,dI(1:3)),L0,L,((n_corr-1)/2+1)*8);
                % ---------------------------------------------------------
                t11t = kt * It;
                
            % WITHOUT CONSTANT TANGENTIAL ACCELERATION  
            else
                t11t = 0;
            end
            
            
            if epsilon_rth && ~flag_eps_r2
            % CONSTANT ACCELERATION IN THE R-T-H FRAME
            % This is formally correct
            % This should be equation [19] - controllare
            t11rth = 3 * kt * cbeta_rth * (calpha_rth *(I(26)-I(25)/(1+P10*s_L0+P20*c_L0))-2*salpha_rth/B*(I(31)-(I(25)*atan((-P10+(P20-1)*tan(L0/2))/B)+B/2*(Ibas(24)*I(25)-corr))));
            % ---------------------------------------------------------
            % Commented by Federico:
            % Stupid old version
            %             t11t=3*kt*cbeta_rth*(calpha_rth*(I(26)-I(25)/(1+P10*s_L0+P20*c_L0))-2*salpha_rth/B*(I(31)+B/2*((I(25)/(1+P10*s_L0+P20*c_L0)-Ibas(24)*I(25)+corr))));
            %                 sqrt(a0^3/muadim)*B*(2*B^2*I(30)+3*P10*I(25))*P11t*0 -...
            %                 sqrt(a0^3/muadim)*B*(2*B^2*I(29)+3*P20*I(25))*P21t*0;
            % ---------------------------------------------------------
            else
                t11rth = 0;
            end
            
        % If initial eccentricity is approximately 0   
        else
            
            if epsilon_t
                % Da dove vengono queste equazioni?
                t11t   = kt*(3/2*(L-L0)^2-4*(1-cos(L-L0)));
                
            else
                t11t = 0;
            end
            
            if epsilon_rth
                t11rth = kt * ( 3 * cbeta_rth * (salpha_rth*((L^2-L0^2)/2-(L-L0)*L0)) -...
                    2*((I(27)*c_L0-I(24)+I(28)*s_L0)*calpha_rth+...
                    4*(-I(27)*s_L0+I(28)*c_L0)*salpha_rth)*cbeta_rth);
            else
                t11rth = 0;
            end
        end
        
        % =================================================================
        % CONSTANT ACCELERATION IN THE INERTIAL REFERENCE FRAME
        % Perche' distingue i seguenti tre casi?
        
        % If initial eccentricity is approximately 0 
        if ~e_flag
            
            if epsilon_In
            t11In = 3 * kt * cbeta0 * ((I(27)-c_L0*I(25)) * cgamma0+...
                (I(28)-s_L0*I(25)) * sgamma0); 
            else
                t11In = 0;
            end
            % -------------------------------------------------------------
            % Commented by Federico:
            % +...
            %                 a0^(3/2)*B/muadim^(1/2)*(- 2*B^2*I(30))*P11t*0 +...
            %                 a0^(3/2)*B/muadim^(1/2)*(- 2*B^2*I(29))*P21t*0;
            % -------------------------------------------------------------
            
        % If eccentricity is not zero but (Omega + omega) is approximately 
        % 90 so that P20 = 0
        elseif ( abs(1/P20) > 1e16)
            if epsilon_In
            t11In = 3 * kt * cbeta0*((I(29)-c_L0/(1+P10*s_L0)*I(25)) * cgamma0-...
                (I(26)-1/(1+P10*s_L0) * I(25))/P10 * sgamma0);
            else
                t11In = 0;
            end
            % -------------------------------------------------------------
            % Commented by Federico:
            % +...
            %                 a0^(3/2)*B/muadim^(1/2)*(- 2*B^2*I(30) - 3*P10*I(25))*P11t*0 +...
            %                 a0^(3/2)*B/muadim^(1/2)*(- 2*B^2*I(29) - 3*P20*I(25))*P21t*0;
            % -------------------------------------------------------------
           
        % If e~=0 and P20~=0
        else
            if epsilon_In
            % Equation [24] - controllata
            t11In = 3 * kt * cbeta0 * ( (-I(26)-P10*I(30)+(1+P10*s_L0)/(1+P10*s_L0+P20*c_L0)*I(25) ) / P20 * cgamma0 + ...
                ( I(30)-s_L0/(1+P10*s_L0+P20*c_L0)*I(25) ) * sgamma0);
            else
                t11In = 0;
            end
            % -------------------------------------------------------------
            % Commented by Federico:
            % +...
            %                 a0^(3/2)*B/muadim^(1/2)*(- 2*B^2*I(30) - 3*P10*I(25))*P11t*0 +...
            %                 a0^(3/2)*B/muadim^(1/2)*(- 2*B^2*I(29) - 3*P20*I(25))*P21t*0;
            % -------------------------------------------------------------
        end
        
        
        % =================================================================
        % eps RTH as 1/r2
        % =================================================================
        if epsilon_rth && flag_eps_r2
            
            n_nodes = max([2 round((L-L0)/(2*pi)*6)]);
            Ls = L0 + (L-L0) * linspace(0, 1, n_nodes+1);
            w  = (L-L0) * [ones(n_nodes-1,1);0.5] / (0.5+n_nodes);
            dfs = dItJ5_dL(Ls(2:end), L0, a0, P10, P20, Q10, Q20, th0, a1rth_r2, P11rth_r2, P21rth_r2);
            t11rth_r2 = B * sqrt(a0 / muadim) * B * quadx_FZ(dfs,w);
            
        else
            t11rth_r2 = 0;
        end
        
        % =================================================================
        % J2 PERTURBATION
        % =================================================================
        if J2
            
            
            % Terms in the time equation:
            % I(39) -> I4c3 =
            % I(40) -> I3c1s3 = 
            % I(35) -> I3c3 = 
            % I(41) -> I2c3s3 = 
            % I(36) -> I2c1s3 = 
            % I(32) -> Icc3
            % I(42) -> I1c3s3
            % I(37) -> I1c2s3
            % I(33) -> Ics3
            % I(29) -> Ic3
            % I(43) -> I4s3
            % I(38) -> I3s3
            % I(34) -> Iss3
            % I(30) -> Is3
            % I(26) -> I13
            
            gt11 = (...
                ((M*(- P10^2*Q10^2 + P10^2*Q20^2 + 4*P10*P20*Q10*Q20 + P20^2*Q10^2 - P20^2*Q20^2))/6) * I(39) +...
                ((4*(P10*Q20 + P20*Q10)*(P10*Q10 - P20*Q20)*M)/3) * I(40) +...
                (((P20*Q10^2 + 2*P10*Q10*Q20 - P20*Q20^2)*M)/2) * I(35) +...
                (-M*(- P10^2*Q10^2 + P10^2*Q20^2 + 4*P10*P20*Q10*Q20 + P20^2*Q10^2 - P20^2*Q20^2)) * I(41) +...
                (-(3*(- P10*Q10^2 + 2*P20*Q10*Q20 + P10*Q20^2)*M)/2) * I(36) +...
                (-((5*P10^2*Q10^4 + 14*P10^2*Q10^2*Q20^2 - 6*P10^2*Q10^2 + 9*P10^2*Q20^4 - 14*P10^2*Q20^2 + P10^2 - 9*P20^2*Q10^4 - 14*P20^2*Q10^2*Q20^2 + 14*P20^2*Q10^2 - 5*P20^2*Q20^4 + 6*P20^2*Q20^2 - P20^2 - 2*Q10^2 + 2*Q20^2))/6) * I(32) +...
                (-(4*(P10*Q20 + P20*Q10)*(P10*Q10 - P20*Q20)*M)/3) * I(42) +...
                (-(3*(P20*Q10^2 + 2*P10*Q10*Q20 - P20*Q20^2)*M)/2) * I(37) +...
                (-(2*(2*P10^2*Q10^3*Q20 + 2*P10^2*Q10*Q20^3 - 4*P10^2*Q10*Q20 - 7*P10*P20*Q10^4 - 14*P10*P20*Q10^2*Q20^2 + 10*P10*P20*Q10^2 - 7*P10*P20*Q20^4 + 10*P10*P20*Q20^2 - P10*P20 + 2*P20^2*Q10^3*Q20 + 2*P20^2*Q10*Q20^3 - 4*P20^2*Q10*Q20 + 2*Q10*Q20))/3) * I(33) +...
                (-((3*P20*Q10^4 - 6*P10*Q10^3*Q20 + 4*P20*Q10^2*Q20^2 - 9*P20*Q10^2 - 6*P10*Q10*Q20^3 - 14*P10*Q10*Q20 + P20*Q20^4 + 5*P20*Q20^2 + 2*P20))/2) * I(29) +...
                ((M*(- P10^2*Q10^2 + P10^2*Q20^2 + 4*P10*P20*Q10*Q20 + P20^2*Q10^2 - P20^2*Q20^2))/6) * I(43) +...
                (((- P10*Q10^2 + 2*P20*Q10*Q20 + P10*Q20^2)*M)/2) * I(38) +...
                (((5*P10^2*Q10^4 + 14*P10^2*Q10^2*Q20^2 - 6*P10^2*Q10^2 + 9*P10^2*Q20^4 - 14*P10^2*Q20^2 + P10^2 - 9*P20^2*Q10^4 - 14*P20^2*Q10^2*Q20^2 + 14*P20^2*Q10^2 - 5*P20^2*Q20^4 + 6*P20^2*Q20^2 - P20^2 - 2*Q10^2 + 2*Q20^2))/6) * I(34) +...
                (-(- P10*Q10^4 + 2*P20*Q10^3*Q20 + 5*P10*Q10^2 + 2*P20*Q10*Q20^3 - 22*P20*Q10*Q20 + P10*Q20^4 - 17*P10*Q20^2)/2) * I(30) +...
                (-((2*P10^2*Q10^4 + 6*P10^2*Q10^2*Q20^2 - 3*P10^2*Q10^2 + 4*P10^2*Q20^4 - 9*P10^2*Q20^2 + P10^2 - 4*P10*P20*Q10^3*Q20 - 4*P10*P20*Q10*Q20^3 + 12*P10*P20*Q10*Q20 + 4*P20^2*Q10^4 + 6*P20^2*Q10^2*Q20^2 - 9*P20^2*Q10^2 + 2*P20^2*Q20^4 - 3*P20^2*Q20^2 + P20^2 + 2*Q10^4 + 4*Q10^2*Q20^2 - 2*Q10^2 + 2*Q20^4 - 14*Q20^2 + 2))/2) * I(26));
            
            gt110 = (-1/12 * (     48 * (P10*Q10 - P10*Q20 - P20*Q10 - P20*Q20) * (P10*Q10 + P10*Q20 + P20*Q10 - P20*Q20) * c_L0^5 +...
                (-96 * (P10*Q10 - P20*Q20) * (P10*Q20 + P20*Q10) * s_L0 - 144 * (P20*Q10^2 + 2*P10*Q10*Q20 - P20*Q20^2)) * c_L0^4 +...
                (-144 * (P10*Q10^2 - 2*P20*Q10*Q20 - P10*Q20^2)*s_L0 + 4*(3*P10^2*Q10^4 - 2*P10^2*Q10^2*Q20^2 - 26*P10^2*Q10^2 - 5*P10^2*Q20^4 + 34*P10^2*Q20^2 - P10^2 - 8*P10*P20*Q10^3*Q20 - 8*P10*P20*Q10*Q20^3 + 80*P10*P20*Q10*Q20 + P20^2*Q10^4 + 2*P20^2*Q10^2*Q20^2 + 6*P20^2*Q10^2 + P20^2*Q20^4 - 14*P20^2*Q20^2 + P20^2 - 28*Q10^2 + 28*Q20^2)) * c_L0^3 +...
                -12*(- P20*Q10^4 + 4*P10*Q10^3*Q20 - 2*P20*Q10^2*Q20^2 - 4*P20*Q10^2 + 4*P10*Q10*Q20^3 - 40*P10*Q10*Q20 - P20*Q20^4 + 12*P20*Q20^2 - P20) * c_L0^2 +...
                12*(P10^2*Q10^4 + 6*P10^2*Q10^2*Q20^2 + 2*P10^2*Q10^2 + 5*P10^2*Q20^4 - 14*P10^2*Q20^2 + P10^2 - 4*P20*P10*Q10*Q20 + Q10^4 + 2*Q10^2*Q20^2 + 2*Q10^2 + Q20^4 - 10*Q20^2 + 1) * c_L0 +...
                8*(4*P10^2*Q10^3*Q20 + 4*P10^2*Q10*Q20^3 - 24*P10^2*Q10*Q20 + P10*P20*Q10^4 - 2*P10*P20*Q10^2*Q20^2 - 10*P10*P20*Q10^2 - 3*P10*P20*Q20^4 + 18*P10*P20*Q20^2 - P10*P20 + 4*P20^2*Q10*Q20 - 28*Q10*Q20) * s_L0^3 +...
                48*(2*P10^2*Q10*Q20 - P10*P20*Q10^4 - P10*P20*Q10^2*Q20^2 + 3*P10*P20*Q10^2 - 2*P10*P20*Q20^2 - P20^2*Q10*Q20 + 4*Q10*Q20) * s_L0 +...
                -6*(P10*Q10^4 - 2*P10*Q10^2*Q20^2 - 10*P10*Q10^2 + 4*P20*Q10*Q20 - 3*P10*Q20^4 + 18*P10*Q20^2 - P10) * sin(2*L0) +...
                -72*Q20*(3*P10*Q10 - P20*Q20)    ) * I(29) +...   & term in I(29) with initial condition
                1/12*((P10^2*Q10^4 + 2*P10^2*Q10^2*Q20^2 + P10^2*Q10^2 + P10^2*Q20^4 - 9*P10^2*Q20^2 + P10^2 - 8*P10*P20*Q10^3*Q20 - 8*P10*P20*Q10*Q20^3 + 20*P10*P20*Q10*Q20 - 5*P20^2*Q10^4 - 2*P20^2*Q10^2*Q20^2 + 19*P20^2*Q10^2 + 3*P20^2*Q20^4 - 11*P20^2*Q20^2 - P20^2 + 28*Q10^2 - 28*Q20^2) * sin(3*L0) +...
                6*(2*P10^2*Q10*Q20 + 3*P10*P20*Q10^4 + 10*P10*P20*Q10^2*Q20^2 - 4*P10*P20*Q10^2 + 7*P10*P20*Q20^4 - 12*P10*P20*Q20^2 + P10*P20 - 4*P20^2*Q10^3*Q20 - 4*P20^2*Q10*Q20^3 + 10*P20^2*Q10*Q20 - 4*Q10*Q20) * c_L0 +...
                -3*(P10^2*Q10^4 + 2*P10^2*Q10^2*Q20^2 - 4*P10^2*Q10^2 + P10^2*Q20^4 - 4*P10^2*Q20^2 + P10^2 - 8*P10*P20*Q10^3*Q20 - 8*P10*P20*Q10*Q20^3 + 24*P10*P20*Q10*Q20 + 15*P20^2*Q10^4 + 22*P20^2*Q10^2*Q20^2 - 32*P20^2*Q10^2 + 7*P20^2*Q20^4 - 8*P20^2*Q20^2 + 3*P20^2 + 4*Q10^4 + 8*Q10^2*Q20^2 - 12*Q10^2 + 4*Q20^4 - 20*Q20^2 + 4) * s_L0 +...
                -6*(3*P20*Q10^4 + 2*P20*Q10^2*Q20^2 - 12*P20*Q10^2 + 8*P10*Q10*Q20 - P20*Q20^4 + 4*P20*Q20^2 + P20) * sin(2*L0) +...
                6*(P10*Q10^4 - 4*P20*Q10^3*Q20 + 2*P10*Q10^2*Q20^2 - 4*P20*Q10*Q20^3 + 16*P20*Q10*Q20 + P10*Q20^4 - 8*P10*Q20^2 + P10) * cos(2*L0) +...
                ((6*P10*P20)*Q10^4 + (-8*P20^2*Q20)*Q10^3 + (4*P10*P20*Q20^2 - 18*P10*P20)*Q10^2 + (10*P10^2*Q20 - 8*P20^2*Q20^3 + 30*P20^2*Q20 + 56*Q20)*Q10 - 2*P10*P20*Q20^4 + 2*P10*P20*Q20^2 + 2*P10*P20) * cos(3*L0) +...
                -18*(P10*Q10^2 - 2*P20*Q10*Q20 - P10*Q20^2) * cos(4*L0) +...
                18*(P20*Q10^2 + 2*P10*Q10*Q20 - P20*Q20^2) * sin(4*L0) +...
                -3*(P10*Q10 - P10*Q20 - P20*Q10 - P20*Q20)*(P10*Q10 + P10*Q20 + P20*Q10 - P20*Q20) * sin(5*L0) +...
                -6*(P10*Q10 - P20*Q20)*(P10*Q20 + P20*Q10) * cos(5*L0) +...
                6*(P10*Q10^4 - 4*P20*Q10^3*Q20 + 2*P10*Q10^2*Q20^2 - 9*P10*Q10^2 - 4*P20*Q10*Q20^3 + 22*P20*Q10*Q20 + P10*Q20^4 + 13*P10*Q20^2 + P10)) * I(30) +...   % term in I(30) with initial condition
                2*((3*Q20^2 - 3*Q10^2)*c_L0^2 + 3*Q10*Q20*sin(2*L0) + (- 2*P20*Q10^2 + 2*P20*Q20^2 - 4*P10*Q10*Q20)*c_L0^3 + (2*P10*Q10^2 - 2*P10*Q20^2 - 4*P20*Q10*Q20)*s_L0^3 + 6*P10*Q10*Q20*c_L0 + 6*P20*Q10*Q20*s_L0) * I(26));
            
            t11J = 3*J2*R^2/(B*G^2*a0)*sqrt(a0/muadim)*(gt11-gt110);

        else
            t11J = 0;
        end
        
        
        % =================================================================
        % J3
        % =================================================================
        if J3
            n_nodes = max([2 round((L-L0)/(2*pi)*6)]);
            Ls = L0 + (L-L0) * linspace(0, 1, n_nodes+1);
            w  = (L-L0) * [ones(n_nodes-1,1);0.5] / (0.5+n_nodes);
            dfs = dItJ5_dL(Ls(2:end), L0, a0, P10, P20, Q10, Q20, th0, a1J3, P11J3, P21J3);
            t11J3 = kJ3 * sqrt(a0 / muadim) * B * quadx_FZ(dfs,w);
            
        else
            t11J3 = 0;
        end
        
        
        % =================================================================
        % J4
        % =================================================================
        if J4
            n_nodes = max([2 round((L-L0)/(2*pi)*6)]);
            Ls = L0 + (L-L0) * linspace(0, 1, n_nodes+1);
            w  = (L-L0) * [ones(n_nodes-1,1);0.5] / (0.5+n_nodes);
            dfs = dItJ5_dL(Ls(2:end), L0, a0, P10, P20, Q10, Q20, th0, a1J4, P11J4, P21J4);
            t11J4 = kJ4 * sqrt(a0 / muadim) * B * quadx_FZ(dfs,w);
            
        else
            t11J4 = 0;
        end
        
        
        
        % =================================================================
        % J5
        % =================================================================
        if J5
            n_nodes = max([2 round((L-L0)/(2*pi)*6)]);
            Ls = L0 + (L-L0) * linspace(0, 1, n_nodes+1);
            w  = (L-L0) * [ones(n_nodes-1,1);0.5] / (0.5+n_nodes);
            dfs = dItJ5_dL(Ls(2:end), L0, a0, P10, P20, Q10, Q20, th0, a1J5, P11J5, P21J5);
            t11J5 = kJ5 * sqrt(a0 / muadim) * B * quadx_FZ(dfs,w);
            
        else
            t11J5 = 0;
        end
        
        
        % =================================================================
        %Sun 
        % =================================================================
        if third_body.flag_Sun
            n_nodes = max([2 round((L-L0)/(2*pi)*6)]);
            Ls = L0 + (L-L0) * linspace(0, 1, n_nodes+1);
            w  = (L-L0) * [ones(n_nodes-1,1);0.5] / (0.5+n_nodes);
            dfs = dItJ5_dL(Ls(2:end), L0, a0, P10, P20, Q10, Q20, th0, a1_3rd_Sun, P11_3rd_Sun, P21_3rd_Sun);
            t11_Sun = eps_3rd_Sun * sqrt(a0 / muadim) * B * quadx_FZ(dfs,w);
            
        else
            t11_Sun = 0;
        end
        
        
        % =================================================================
        % Moon 
        % =================================================================
        if third_body.flag_Moon
            n_nodes = max([2 round((L-L0)/(2*pi)*6)]);
            Ls = L0 + (L-L0) * linspace(0, 1, n_nodes+1);
            w  = (L-L0) * [ones(n_nodes-1,1);0.5] / (0.5+n_nodes);
            dfs = dItJ5_dL(Ls(2:end), L0, a0, P10, P20, Q10, Q20, th0, a1_3rd_Moon, P11_3rd_Moon, P21_3rd_Moon);
            t11_Moon = eps_3rd_Moon * sqrt(a0 / muadim) * B * quadx_FZ(dfs,w);
            
        else
            t11_Moon = 0;
        end
        
        % =================================================================
        % DRAG
        % =================================================================
        if drag.CD
            t11Drag = It1_Drag;
        else
            t11Drag = 0;
        end
        


        % =================================================================
        % Time
        % =================================================================
        t(i) = t00 + ...
            epsilon_t * t11t + ...
            epsilon_rth * t11rth + ...
            epsilon_rth * t11rth_r2 + ...
            epsilon_In * t11In + ...
            t11J + ...
            t11J3 + ...
            t11J4 + ...
            t11J5 + ...
            t11_Sun + ...
            t11_Moon + ...
            drag.eps * t11Drag;

    end
%     keyboard
%     if t(i) < -1e-15
%         warning('Negative time in analytical propagation')
%         keyboard
%     end
    

    %% Calculate new values for the Equinoctial elements
%     keyboard
    Equin(:,i) =              [a0;    P10;    P20;    Q10;    Q20;    L + (Equin0(6)-L0)] + ...        % Unperturbed equinoctial elements
                  epsilon_t  *[a1t;   P11t;   P21t;   Q11t;   Q21t;   0] + ...                        % Constant tangential acceleration
                  epsilon_rth*[a1rth; P11rth; P21rth; Q11rth; Q21rth; 0] + ...                        % Constant acceleration in the rth frame
                  epsilon_rth*[a1rth_r2; P11rth_r2; P21rth_r2; Q11rth_r2; Q21rth_r2; 0] + ...         % Constant acc in the rth frame with eps as 1/r2
                  epsilon_In *[a1In;  P11In;  P21In;  Q11In;  Q21In;  0] + ...                         % Constant intertial acceleration
                              [a1J;   P11J;   P21J;   Q11J;   Q21J;   0] + ...                         % J2 acceleration
                              [a1J3;  P11J3;  P21J3;  Q11J3;  Q21J3;  0] + ...                         % J3 acceleration     
                              [a1J4;  P11J4;  P21J4;  Q11J4;  Q21J4;  0] + ...                         % J4 acceleration 
                              [a1J5;  P11J5;  P21J5;  Q11J5;  Q21J5;  0] + ...                         % J5 acceleration 
                  eps_3rd_Sun  *  [a1_3rd_Sun; P11_3rd_Sun; P21_3rd_Sun; Q11_3rd_Sun; Q21_3rd_Sun; 0] + ...                       % third body
                  eps_3rd_Moon  *  [a1_3rd_Moon; P11_3rd_Moon; P21_3rd_Moon; Q11_3rd_Moon; Q21_3rd_Moon; 0] + ...                       % third body 
                  drag.eps *  [a1Drag; P11Drag; P21Drag; Q11Drag; Q21Drag; 0];                        % Drag acceleration     

%  keyboard
% 
else
    % The difference between the final and initial true longitude is not
    % big enough. Return the initial equinoctial elements and the initial
    % time
    Equin = Equin0;
    t = 0;
end


end

% if Equin(1) > 10
%     keyboard
% end
% P10
% P11Drag
% drag.CD * drag.A_m * a0^2 * B^2 / 2

return

% 