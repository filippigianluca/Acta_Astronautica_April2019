function [Equin,t] = AnEquin_all_forward_tang_om_change_in_thrust_multi_m2(L, Equin0, ...
                                                                          epsilon, beta_t, u_ratio, ...
                                                                          epsilon_rth, alpha_rth, beta_rth, ...
                                                                          epsilon_Sun, alpha_Sun, beta_Sun, ...
                                                                          muadim, geopotential, third_body, R, drag, ...
                                                                          Earth_flat,ill_flag, n_FPET, ...
                                                                          flag_eps_r2)

                                                                      
% Input:

% Output:


% Author: Federico Zuiani
% Comments and code tyding: Marilena Di Carlo

if nargin < 20
    flag_eps_r2 = 0;
end

dL = ( L(2)-L(1) ) / n_FPET;
L_t = L(1) + dL;

% Keplerian elements corresponding to initial equinoctial elements
kep0 = eq2kep(Equin0);

% Rotation matrix (rotation around x axis, angle equal to inclination)
Rb00 =  [1       0             0; ...
        0 cos(kep0(3)) -sin(kep0(3));...
        0 sin(kep0(3))  cos(kep0(3))];
  
% Rotation matrix (rotation around z axis, angle equal to omega + f, where
% omega is the perigee argument and f is the true anomaly)
Rc00 = [cos(Equin0(6)-kep0(4)) -sin(Equin0(6)-kep0(4)) 0;...
        sin(Equin0(6)-kep0(4))  cos(Equin0(6)-kep0(4)) 0;...
                  0                      0             1];

Eq_t = Equin0;
tt = 0;

for i = 1 : n_FPET
    
    % Square of eccentricity
    e_2 = Eq_t(2)^2 + Eq_t(3)^2;
    
    % e cos(theta) where theta = L - (Omega + omega) is the true anomaly
    e_c_theta = Eq_t(2) * sin(Eq_t(6)) + Eq_t(3) * cos(Eq_t(6));
    e_s_theta = Eq_t(3) * sin(Eq_t(6)) - Eq_t(2) * cos(Eq_t(6));
    
    D = sqrt(e_2 + 1. + 2*(e_c_theta));
    
    % u_ratio is 0 for energy encreasig law, -1 for perigee argument
    % decreasing law and 1 for perigee argument increasing law
    %
    % If we are considering a perigee argument decreasing or increasing law
    if abs(u_ratio) > 0
        
        % ????
        alfa_om = atan2(-e_s_theta, e_c_theta * (1+e_c_theta) / (2+e_c_theta) ) + pi;
        
        if u_ratio<0
            alfa_om = alfa_om-pi;
            u_ratio = -u_ratio;
        end
        
        %% WIP
        % -----------------------------------------------------------------
        % Commented by Federico
        %     e=sqrt(e_2);
        %     c_theta=e_c_theta/e;
        %     s_theta=e_s_theta/e;
        %     k_tang_max=(1+e);
        %     k_tang_min=(1-e);
        % -----------------------------------------------------------------
        k_tang = 1;
        % -----------------------------------------------------------------
        % Commented by Federico
        %0.5+0.5*(D-k_tang_min)/(k_tang_max-k_tang_min);%1;%
        %     O=e*(9-9*e_2+sqrt(3)*sqrt(27-54*e_2+27*e_2^2+4*e_2^3))^(1/3);
        %     th_max=acos(-(1/e) - ((2/3)^(1/3)*e_2)/O+O/(2^(1/3)*3^(2/3)*e_2));
        %     k_om_max=abs((-cos(th_max)-((2+e*cos(th_max))^2*sin(th_max)*tan(th_max))/(1+e*cos(th_max))^2)/sqrt(1+(e*sin(th_max)+2*tan(th_max))^2/(1+e*cos(th_max))^2));
        %     k_om_min=1;
        %     DD=abs((-c_theta-((2+e_c_theta)^2*s_theta^2/c_theta)/(1+e_c_theta)^2)/sqrt(1+(e_s_theta+2*s_theta/c_theta)^2/(1+e_c_theta)^2));
        % -----------------------------------------------------------------
        k_om = 1;
        % -----------------------------------------------------------------
        % Commented by Federico
        %0.5+0.5*(DD-k_om_min)/(k_om_max-k_om_min);%1;%
        %     if u_ratio>=0.5
        %         plot(tt,[k_tang;k_om;k_tang/k_om],'.')%
        %     end
        % -----------------------------------------------------------------
        %%
        
        w_tang = (1-u_ratio) * k_tang;
        w_om   = abs(u_ratio) * k_om;
        
        if epsilon_Sun
            
            if i==1
                u_Sun = [cos(beta_In) * cos(alpha_In); ...
                         cos(beta_In) * sin(alpha_In); ...
                         sin(beta_In)];
            else
                kep_t = eq2kep(Eq_t);
                
                Rb=[1 0 0;...
                    0 cos(kep_t(3)) sin(kep_t(3));...
                    0 -sin(kep_t(3)) cos(kep_t(3))];
                
                Rc=[cos(Eq_t(6)-kep_t(4)) sin(Eq_t(6)-kep_t(4)) 0;...
                    -sin(Eq_t(6)-kep_t(4)) cos(Eq_t(6)-kep_t(4)) 0;...
                    0 0 1];
                
                R1=[cos(kep_t(4)-kep0(4)) sin(kep_t(4)-kep0(4)) 0;...
                    -sin(kep_t(4)-kep0(4)) cos(kep_t(4)-kep0(4)) 0;...
                    0 0 1];
                
                u_Sun=(Rc*Rb*R1*Rb00*Rc00)*[cos(beta_Sun)*cos(alpha_Sun);...
                                            cos(beta_Sun)*sin(alpha_Sun);...
                                            sin(beta_Sun)];
            end
        else
            u_Sun=[0;0;0];
        end
        
        u_r = w_tang * cos(beta_t) * e_s_theta/D+w_om*cos(alfa_om);
        u_t = w_tang * cos(beta_t) * (1+e_c_theta)/D+w_om*sin(alfa_om);
        u_h = w_tang * sin(beta_t);
        u = sqrt(u_r^2+u_t^2+u_h^2);
        
        w_tang = w_tang/u;
        w_om   = w_om/u;
        
        u_In = [epsilon*w_om*cos(alfa_om)+epsilon_Sun*u_Sun(1);...
            epsilon*w_om*sin(alfa_om)+epsilon_Sun*u_Sun(2);...
            epsilon_Sun*u_Sun(3)];
        
        epsilon_In=sqrt(u_In(1)^2+u_In(2)^2+u_In(3)^2);
        beta_In=asin(u_In(3)/epsilon_In);
        alpha_In=atan2(u_In(2),u_In(1));
        
        
    % Is u_ratio is equal to zero, meaning an energy increasing law
    else
        
        w_tang = 1;
        
        if epsilon_Sun
            epsilon_In = epsilon_Sun;
            
            if i == 1
                alpha_In = alpha_Sun;
                beta_In  = beta_Sun;
            else
                
                % I think that here an update of the Sun angles is applied
                kep_t = eq2kep(Eq_t);
                
                % Rotation matrix - rotation around x axis of an angle
                % equal to the inclination in counterclockwise direction
                Rb = [1 0 0; ...
                      0 cos(kep_t(3)) sin(kep_t(3));...
                      0 -sin(kep_t(3)) cos(kep_t(3))];
                 
                % Rotation matrix -   
                Rc = [cos(Eq_t(6)-kep_t(4)) sin(Eq_t(6)-kep_t(4)) 0; ...
                     -sin(Eq_t(6)-kep_t(4)) cos(Eq_t(6)-kep_t(4)) 0; ...
                     0 0 1];
                 
                 
                R1 = [cos(kep_t(4)-kep0(4)) sin(kep_t(4)-kep0(4)) 0;-sin(kep_t(4)-kep0(4)) cos(kep_t(4)-kep0(4)) 0;0 0 1];
                u_Sun = (Rc*Rb*R1*Rb00*Rc00)*[cos(beta_Sun)*cos(alpha_Sun);cos(beta_Sun)*sin(alpha_Sun);sin(beta_Sun)];
                beta_In = asin(u_Sun(3));
                alpha_In = atan2(u_Sun(2),u_Sun(1));
            end
        else
            beta_In=0;
            alpha_In=0;
            epsilon_In=0;
        end
        
    end
    
%     [Equin,t_] = AnEquin_all_forward_tang_NoDrag(L_t, Eq_t, epsilon*w_tang, beta_t, epsilon_rth, alpha_rth, beta_rth, epsilon_In, alpha_In, beta_In, muadim, J2, R, ill_flag);
    [Equin,t_] = AnEquin_all_forward_tang_m(L_t, Eq_t, epsilon*w_tang, beta_t, epsilon_rth, alpha_rth, beta_rth, epsilon_In, alpha_In, beta_In, ...
        muadim, geopotential, third_body, R, drag, Earth_flat, ill_flag, flag_eps_r2);

    % Redefine Eq_t as final point of propagation
    Eq_t = Equin;
    
    tt=tt+t_;
    L_t=L_t+dL;
end

t=tt;

return