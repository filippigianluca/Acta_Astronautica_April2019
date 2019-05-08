function [I,dI] = Integrals_tang_m(L0, L, P10, P20, e_flag)

% Calculates analytical integrals for the perturbative expansion of the
% constant tangential acceleration case. 
%
% INPUTS:
% L0: longitude of the initial state.
% L:  longitude of the current state, can be a vector.
% P10, P20: values at x(1) for Equinoctial elements P1 and P2.
% e_flag: 
%
% OUTPUTS:
% I: Integrals Ia, Ip1 and Ip2 of the Constant Tangential Acceleration case
        % The output of I is, for each value of L, is:
        % I = [Ia Ip1 Ip2] where Ip1 and Ip2 are not the one in the paper
        % but:
        % IP1 = s_Omom*IPa + c_Omom*IPb;
        % IP2 = c_Omom*IPa-s_Omom*IPb;
        % and IPa and IPb are what in the paper are called Ip1 and Ip2
% dI:

% REFERENCE:
% Zuiani, "Extended Analytical Formulas for the Perturbed Keplerian Motion 
% Under a Constant Control Acceleration"

% Author: Federico Zuiani
% Comments and code tidying: Marilena Di Carlo



if any(abs(L-L0)>1e-6)
    
    [a,b]=size(L);

    if a > b
        L=L.';
    end

    e0 = sqrt(P10^2 + P20^2);
    errtol = 1e-8;
    
    % If initial eccentricity is not zero
    if e_flag
        % cos(Omega+omega)
        c_Omom = P20/e0;
        
        % sin(Omega+omega)
        s_Omom=P10/e0;
        
        % -----------------------------------------------------------------
        % Commented by Federico:
        %         Omom=atan2(s_Omom,c_Omom);
        % -----------------------------------------------------------------
        
        c_L0 = cos(L0);
        s_L0 = sin(L0);
        c_L  = cos(L);
        s_L  = sin(L);
        
        % e*cos(theta0) where L = (theta + Omega + omega)
        e0_c_th0 = P10*s_L0+P20*c_L0;
        
        % e*sin(theta0) where L = (theta + Omega + omega)
        e0_s_th0 = P20*s_L0-P10*c_L0;
        
        % theta0
        th0 = atan2(e0_s_th0,e0_c_th0);
        
        % -----------------------------------------------------------------
        % Commented by Federico:
        %         th0=L0-Omom;
        %         e0_c_th0=e0*cos(th0);
        %         e0_s_th0=e0*sin(th0);
        % -----------------------------------------------------------------
        
        % e*cos(theta) where L = (theta + Omega + omega)
        e0_c_th = P10*s_L+P20*c_L;
        
        % e*sin(theta) where L = (theta + Omega + omega)
        e0_s_th = P20*s_L-P10*c_L;
        
        % theta
        th=atan2(e0_s_th,e0_c_th);
        
        % -----------------------------------------------------------------
        % Commented by Federico:
        %         th=L-Omom;
        %         e0_c_th=e0*cos(th);
        %         e0_s_th=e0*sin(th);
        % -----------------------------------------------------------------
        
        % Attention! The following H and H0 have not the same meaning as in
        % the others Integrals functions
        H0 = sqrt(1+e0^2+2*e0_c_th0);
        H  = sqrt(1+e0^2+2*e0_c_th);
        
        K0 = 1+e0_c_th0;
        K = 1+e0_c_th;
        
        % Elliptic integrals
        ElliEtmp = lellipe3([th0 th pi]/2,4*e0/(1+e0)^2,errtol);
        ElliFtmp = lellipf3([th0 th pi]/2,4*e0/(1+e0)^2,errtol);
        
        % -----------------------------------------------------------------
        % Commented by Federico:
        %         tic
        %         toc
        %         tic
        %         ElliFtmp=lellipf3([th0 th]/2,4*e0/(1+e0)^2,errtol);
        %         toc
        %         tic
        %         ElliFtmp=lellipf3([th0 th]/2,4*e0/(1+e0)^2,errtol);
        %         toc
        %         tic
        %         ElliFtmp=lellipf3([th0 th]/2,4*e0/(1+e0)^2,errtol);
        %         toc
        %         tic
        %         ElliFtmp=lellipf3([th0 th]/2,4*e0/(1+e0)^2,errtol);
        %         toc
        %         tic
        %         ElliFtmp=lellipf3([th0 th]/2,4*e0/(1+e0)^2,errtol);
        %         toc
        
        %         ElliE=0*L;
        %         ElliF=ElliE;
        %         for i=1:length(L)
        %             ElliE(i)=lellipe2(th(i)/2,4*e0/(1+e0)^2,errtol)-ElliE0;
        %             ElliF(i)=lellipf2(th(i)/2,4*e0/(1+e0)^2,errtol)-ElliF0;
        %         end
        % -----------------------------------------------------------------
        
        ElliE = ElliEtmp(2:end-1) - ElliEtmp(1);
        ElliF = ElliFtmp(2:end-1) - ElliFtmp(1);
        
        % Useful precomputation for Equation [15] and [16] in the Appendix:
        
        % Arcotangent term in equation [16] in Appendix. (ok)
        % Federico uses the following relationship for difference of arcotangent:
        % atan(x1) - atan(x2) = atan((x1 - x2) / (1 +x1 x2))
        Atan = atan( sqrt(1-e0^2)*(H-H0)./((1-e0^2) + H*H0) );
        
        % -----------------------------------------------------------------
        % Commented by Federico:
        %atan(H/sqrt(1-e0^2))-atan(H0/sqrt(1-e0^2));%
        % -----------------------------------------------------------------
        
        % This is the logarithmic term in Equation [15], considering both
        % theta and theta0 (controllato)
        Atanh = 0.5 * log( (H0 - e0_s_th0) * (H + e0_s_th) ./ ( (H-e0_s_th) * (H0+e0_s_th0) ) );
        % -----------------------------------------------------------------
        % Commented by Federico:
        %atanh(e0_s_th/H)-atanh(e0_s_th0/H0);%
        % -----------------------------------------------------------------
        
        % Ia (ok)
        Ia = ( ElliE/(1-e0) + ElliF/(1+e0) - (H.*e0_s_th./K - H0*e0_s_th0/K0) / (1-e0^2));

        
        % -----------------------------------------------------------------
        % Commented by Federico:
        %         IPa=(ElliE/(e0*(1-e0))-(1+e0^2)*ElliF/(e0*(1+e0)^2*(1-e0))+Atanh/(e0*(1-e0)*(1+e0))-(H.*e0_s_th./K-H0*e0_s_th0/K0)/(1+e0))/(e0*(1-e0));
        %         IPa=ElliE/(e0*(1-e0))-(1+e0^2)*ElliF/(e0*(1-e0)*(1+e0)^2)+1/(e0*(1-e0)*(1+e0))*(atanh(e0*sin(th)./sqrt(1+e0^2+2*e0*cos(th)))-atanh(e0*sin(th0)./sqrt(1+e0^2+2*e0*cos(th0))))-(sqrt(1+e0^2+2*e0*cos(th)).*(1+cos(th)).*tan(th/2)./(1+e0*cos(th))-sqrt(1+e0^2+2*e0*cos(th0)).*(1+cos(th0)).*tan(th0/2)./(1+e0*cos(th0)))/((1+e0)*(1-e0));
        % -----------------------------------------------------------------
       
        % Ip1
        % The following is the integral in Equation [15] in the Appendix
        % (ok)
        IPa=( ElliE - (1+e0^2)/(1+e0)^2 * ElliF + Atanh/(1+e0) - (H.*e0_s_th./K - H0*e0_s_th0/K0) / (1+e0)) / (e0*(1-e0));
        
        % The following is the integral in Equation [16] in the Appendix
        % (ok)
        IPb = -1/e0 * ( 2*Atan/(1-e0^2)^(3/2) + (H./K-H0/K0)/(1-e0^2));
        
        IP1 = s_Omom*IPa + c_Omom*IPb;
        
        % Ip2
        IP2 = c_Omom*IPa-s_Omom*IPb;
        
        % =================================================================    
        % Discontinuities at pi (same as in Limits4_m.m for the other
        % integrals)
        % dI will be summed to the analytic integrals computed aboce
        
        DIPa=( ElliEtmp(end) - (1+e0^2) / (1+e0)^2 * ElliFtmp(end) ) / (e0*(1-e0));

        % =========================  First term
        % For Ia the value of the integral at pi is given only by the terms
        % in EllipticE and EllipticF (that is why above the elliptic
        % integrals are computed also for pi!). The cofficient 2 comes from
        % the fact that the left and right limits are equal (otherwise dI
        % would be zero)
        % For the second and third term remember that Ip2 as reported in
        % Equation 16 has Ip2[pi] = Ip2[-pi] (that is why only the term
        % DIpa appears)
        % ========================= Second term
        % Remember that in this file it is not like in the paper; the
        % second output is sin(Omega+omega) Ip1 + cos(Omega+omega)Ip2 as in
        % Equation 28 and Ip1 and Ip2 are equations 15 and 16 in Appendix
        % ========================= Third term
        % Remember that in this file it is not like in the paper; the
        % second output is cos(Omega+omega) Ip1 - sin(Omega+omega)Ip2 as in
        % Equation 28 and Ip1 and Ip2 are equations 15 and 16 in Appendix
        dI=[ 2 * (ElliEtmp(end)/(1-e0) + ElliFtmp(end)/(1+e0)); ...
             2 * P10 * DIPa / e0; ...
             2 * P20 * DIPa / e0];
        
    % If e0 = 0, then  
    else
        
        % If e0 = 0, then the integrand of Ia is 1
        Ia  = L - L0;
        
        % If e0=0, then:
        % Ipa = sin(theta) - sin(theta0)
        % Ipb = cos(theta0) - cos(theta)
        % 
        % and therefore, considering Equation [28b] we obtain:
        % 
        % Ip1 = sin(Omega+omega) sin(theta) - sin(Omega+omega) sin(theta0)
        % + cos(Omega+omega) cos(theta0) - cos(Omega+omega) cos(theta) = 
        % = - [ cos(Omega+omega) cos(theta) - sin(Omega+omega) sin(theta)]
        % + [cos(Omega+omega) cos(theta0) - sin(Omega+omega) sin(theta0)] =
        % = - cos L + cos L0 
        % 
        % Ip2 = cos(Omega+omega) sin(theta) - cos(Omega+omega) sin(theta0)
        % - sin(Omega+omega) cos(theta0) + sin(Omega+omega) cos(theta) = 
        % = sinL - sin L0
        
        IP1 = -(cos(L)-cos(L0));
        IP2 = sin(L)-sin(L0);
        
        dI  = [0;0;0];
        
    end
    
    I = [Ia; IP1; IP2];
else
    
    I = zeros(3,length(L));
    dI = [0;0;0];
    
end

return