function [dx_dt, t_f, x_f] = integrand_analytic_propagator_new(t, x, control, input, constants)                                 
   
% Function for the integration of the averaged orbital elements. 

% IMPORTANT!!!!!!
% when integrand_analytic_propagator uses perigee and apogee arcs, the
% thrust on those two arcs has azimuth angle equal to zero. It is not
% possible to set the azimuth angle on those two thrust arcs.
% This function allows to set the value of the thrust azimuth on the
% perigee and apogee thrust arcs (in a similar way to beta). 
% The previous function was not modified but a new function was created to
% allow to define alpha on the perigee and apogee thrust arc. 

% Also, different method to compute the shift of the perigee thrust arc wrt
% to the perigee than Federico
%%%%%%%%%%%%%%%%%%%%%%%%

% =========================================================================
% INPUT
% =========================================================================

% t   ->
% x   -> Vector composed by equincotial elements(row 1-6), mass(row 7), 
%        time spent below the radiation belt (row 8) and time in eclipse
%        condition (row 9)
% control   -> structure with THRUST VECTOR CONTROL PARAMETERS
% input     -> structure with INPUT
% constants -> structure with CONSTANTS

% ---------------THRUST VECTOR CONTROL PARAMETERS -------------------------
% dL_a -> Semi-amplitude of the APOCENTRE thrusting arc, in rad between -pi and pi
% dL_p -> Semi-amplitude of the PERICENTRE thrusting arc, in rad between -pi and pi
% k_a  -> k_a = 1, apocenter is raised
% k_p  -> k_p = 1, pericenter is raised
% alfa -> alfa is defined equal to pi/2 above, don't know why...
% beta_a -> Angle between the orbit plane and the out-of-plane thrust 
%           component, in rad, between -pi/2 and pi/2.
% beta_p -> Same as beta_a, but for the pericenter thrusting
% eta    -> Offset angle of the midpoint of the pericentre thrusting arc
%           w.r.t. the pericentre itself, in rad, positive counter-clockwise. 
% csi    -> Ratio between the amplitudes of the post-apocentre and 
%           post-pericentre coasting arcs. It is comprised between 0 and 1.
% u_ratio -> Defines the ratio between the energy-increasing and
%            small-omega-changing, locally optimal, control laws. 
%            u_ratio=0 means a purely energy-increasing law. 
%            u_ratio=1 means a purely small-omega-increasing law. 
%            u_ratio=-1 means a purely small-omega-decreasing law.
% -------------------------------------------------------------------------

% -------------------  CONSTANTS ------------------------------------------
% mu     -> gravitational parameter, adimensional
% J2, J3, J4, J5 
% R_Earth-> Earth radius

% -----------------INPUT: LOW THRUST ENGINE PARAMETERS---------------------
% T_adim -> Maximum engine thrust 
% m_rate -> Specific mass consumption = 1/(I_sp*g_0)

% -------------------INPUT: OTHER INPUTS ----------------------------------
% T_Sun_adim -> Magnitude of the SRP thrust
% e_target -> Targeted final value of the eccentricity
% r_belt   -> Radius of the radiation belt
% e_flag   -> Eclipse flag:   0 - no eclipse are considered
%                             1 - eclipses are accounted for

% ---------------------- TIMES --------------------------------------------
% t_0 -> Reference initial time + Waiting time
% t_om_change -> Reference initial time + Waiting time + Time elapsed, from 
%                the end the commissioning phase, from which a purely 
%                small-omega-changing control law is used.

% Author: Federico Zuiani, Marilena Di Carlo
% marilena.di-carlo@strath.ac.uk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants

muadim = constants.mu;
R      = constants.R_Earth;
g0     = constants.g0;


%% Initialization 

% Control variables
dL_a    = control.dL_a;
dL_p    = control.dL_p;
k_a     = control.k_a;
k_p     = control.k_p;
alfa_a  = control.alfa_a;
alfa_p  = control.alfa_p;
beta_a  = control.beta_a;
beta_p  = control.beta_p;
eta     = control.eta;
csi     = control.csi;
u_ratio = control.u_ratio;
ts      = control.ts;

% Input
drag         = input.drag;
thrust       = input.thrust;
T_adim       = input.T_adim;
m_rate       = input.m_rate;
T_Sun_adim   = input.T_Sun_adim;
r_belt       = input.r_belt;
t_0          = input.t_0;
options      = input.options_SOL;
ecl_flag     = input.ecl_flag;
Earth_flat   = input.Earth_flattening;

J2 = input.geopotential.J2;

geopotential = input.geopotential;


% Equinoctial elements
Eq_a_o = x(1:6);

% Mass
m = x(7);

if any(isnan(x))
    dx_dt = 0 * x;
    return
end


% =========================================================================
% Interpolation
% =========================================================================
% Interpolation of dL_a over time going from 0 to maximum time ttt (ttt and
% dL_a have the same number of elements)
% interp1q_c is a mex file
% tt is a vector of times going from 0 to T_max_adim and with n_x elements

% interp1q_c works as follow: a series of apogee true longitude dL_a are defined
% at the times tt. interp1q_c computes the value of dL_a at times ttt, that
% is dL_a_i.
% The same for all the other variable
dL_a_i   = interp1q_c(ts, dL_a, abs(t));
dL_p_i   = interp1q_c(ts, dL_p, abs(t));
beta_a_i = interp1q_c(ts,beta_a,abs(t));
beta_p_i = interp1q_c(ts,beta_p,abs(t));
alfa_a_i = interp1q_c(ts,alfa_a,abs(t));
alfa_p_i = interp1q_c(ts,alfa_p,abs(t));
eta_i    = interp1q_c(ts,eta,abs(t));
csi_i    = interp1q_c(ts,csi,abs(t));
u_ratio_i = interp1q_c(ts,u_ratio,abs(t));


% ---------------- DECOMMENTARE -------------------------------------------
% % Non ho idea di cosa fa qui - per ora lo commento perche' relativo alla
% % parte di controllo con thrust
% if (t_0 + ttt * t/86400  >=  t_om_change )  &&  ( t_0 + ttt * t/86400 < t_om_change + abs(Dt_om_change) )
%     if Dt_om_change<0
%         u_ratio_i=-1;
%     else
%         u_ratio_i=1;
%     end
% end
% -------------------------------------------------------------------------

% If apogee arcs are lower than 0, switch k_a (no more apogee raising but
% apogee decrease), switch the values of dL_a_i and also beta_a_i
% (why?!?!?)
if dL_a_i < 0
    k_a      = -k_a;
    dL_a_i   = -dL_a_i;
    beta_a_i = -beta_a_i;
end

% If perigee arcs are lower than 0, switch k_p (no more perigee raising but
% perigee decrease), switch the values of dL_p_i and also beta_p_i
% (why?!?!?)
if dL_p_i < 0
    k_p      = -k_p;
    dL_p_i   = -dL_p_i;
    beta_p_i = -beta_p_i;
end

% Is sum of apogee and perigee arc is greater than pi, reduce both dL_a_i
% and dL_p_i in such a way that their sum is pi at maximum
if dL_a_i + dL_p_i > pi
    dL     = dL_a_i + dL_p_i;
    dL_a_i = pi * dL_a_i/dL;
    dL_p_i = pi*dL_p_i/dL;
end

% Non so come questo possa succedere dopo il precedente ciclo if, comunque
if dL_a_i > pi
    dL_a_i = pi;
    dL_p_i = 0;
elseif dL_p_i > pi
    dL_p_i = pi;
    dL_a_i = 0;
end

% Remember that dL_p_i and dL_a_i are semiamplitude. Therefore the
% following lines computes the remaining arc over one period:
% This is DeltaL_p + DeltaL_a
% The 2*pi angle is indeed thus subdivided: 2*dL_p (pericenter thrust angle),
% DeltaL_p (post-pericenter coasting angle), 2*dL_a (apocenter thrust
% angle), DeltaL_a (post-apocenter coasting angle) and they are such that:
% 2*dL_p + DeltaL_p + 2*dL_a + DeltaL_a = 2*pi
% dL_c = DeltaL_p + DeltaL_a
dL_c = 2 * (pi - (dL_p_i+dL_a_i) );

% This SHOULD be the post-apocenter thrusting arc amplitude, if csi is 
% csi = DeltaL_a / (DeltaL_a + DeltaL_p)
dL_c1_i=csi_i*dL_c;

% -------------------------------------------------------------------------
% Commented by Federico:
% dL_c2_i=dL_c-dL_c1_i;
% -------------------------------------------------------------------------

% Cosa fa qui?
ill_flag = 0;
if abs(u_ratio_i)>0.20
    dL_subarc=2*pi/4;
else
    dL_subarc=2*pi/2;
end

t_ecl=0;
th_ecl_med=0;
t_ecl_switch=1;
adj_flag=0;
DT=0;

% Choosing correct apses point ("a" does not necessarily mean apocenter nor
% "p" pericenter, since they could swap positions if the zero eccentricty
% threshold is crossed)

% Initial keplerian elements
kep_in = eq2kep(x(1:6));


if abs(kep_in(6)-pi)<min([kep_in(6) 2*pi-kep_in(6)])
    th_in=kep_in(4)+kep_in(5)+pi;
else
    th_in=kep_in(4)+kep_in(5);
end



%% Sun position for third body

if input.flag_3rd_Sun
    
    % Earth Sun position in the ecliptic plane
    rr = -EphSS(3, t_0 + t*constants.TU/86400);
    
    % Ecliptic inclination
    incl_ecl = astro_constants(8);
    
    % Earth Sun position in the ECI reference frame
    r_E_S = [rr(1); cos(incl_ecl)*rr(2); sin(incl_ecl)*rr(2)] ./ constants.DU;
    
    %
    third_body.flag_Sun         = 1;
    third_body.R_3rd_Sun_vector = r_E_S;
    third_body.mu_3rd_Sun       = constants.mu_Sun;
    
else
    third_body.flag_Sun = 0;
    
end


%% Moon position for third body

if input.flag_3rd_Moon
    
    third_body.flag_Moon = 1;
    third_body.R_3rd_Moon_vector = EphSS(11, t_0 + t*constants.TU/86400) / constants.DU;
    third_body.mu_3rd_Moon = constants.mu_Moon;
else
    third_body.flag_Moon = 0;
    
end

%% Coasting after apocenter with J2 and SRP (if not in eclipse) perturbation

% keyboard

% True longitude of the initial equinoctial elements point after the
% apocenter thrust arc. Start from pericenter (Eq_a_o(6) correspond to
% true anomaly equal to zero), add shift from pericenter for the center of
% the pericenter thrust arc, then subtract half of the pericenter angle and
% then again substract post-apocenter thrusting arc. You will end up at the
% end of the apocenter thrusting arc
Eq_a_o(6) = Eq_a_o(6) + eta_i - dL_p_i - dL_c1_i;
Eq_a_o(6) = 0 + eta_i + atan(Eq_a_o(2)/Eq_a_o(3)) - dL_p_i;
Eq_a_o(6) = mod(Eq_a_o(6),2*pi);

Eq_i = Eq_a_o;

% keyboard
% Adjust L_p_i such that it is increasing w.r.ttt. Eq_a_o
% L_p_i is the true longitude corresponding to the beginning of the perigee 
% thrust arc
L_p_i = Eq_a_o(6) + dL_c1_i;


% -------------------------------------------------------------------------
% Commented by Federico:
% while (L_p_i-Eq_a_o(6))>=2*pi
%     L_p_i=L_p_i-2*pi;
% end
% while (L_p_i-Eq_a_o(6))<0
%     L_p_i=L_p_i+2*pi;
% end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Check shadow position and eventually adjust beginning of thrusting arc
% -------------------------------------------------------------------------
% If ecl_flag is zero we are not taking into account the eclipses
if ecl_flag == 0
    
    L_i = NaN;
    L_o = NaN;
    alpha_Sun = 0;
    beta_Sun  = 0;
    phi_Sun   = 0;
    
else
    
    % L_i = true longitude eclipse entry point
    % L_o = true longitude eclipse exit point
    % kep: keplerian parameters in the J2000, Earth-centered equatorial frame.
    % dt_ecl: eclipse duration.
    % kep2: keplerian parameters in an Earth-centered, Ecliptic reference frame
    %       in which the x direction is aligned with the Earth-Sun vector.
    % alpha_Sun: azimuth of the Sun-Earth direction (i.e. of the solar
    %            radiation pressure acceleration vector) in the r-t-h reference frame.
    % beta_Sun: elevation of the Sun-Earth direction (i.e. of the solar
    %           radiation pressure acceleration vector) in the r-t-h reference frame.
    % phi_Sun: 

    [L_i, L_o, kep, dt_ecl, kep2, alpha_Sun, beta_Sun, phi_Sun] = sol_Eclipse_Eq(Eq_i, t_0 + t*constants.TU/86400, constants.mu, constants.R_Earth);
%     [L_i, L_o, kep, dt_ecl, kep2, alpha_Sun, beta_Sun, phi_Sun] = sol_Eclipse_Eq(Eq_i, t_0 + ttt, muadim, R);

end

% -------------------------------------------------------------------------
% Commented by Federico:
% [Eq_a_o(6) L_p_i L_i L_o]
% -------------------------------------------------------------------------

% Non capisco cosa fa se non so cos'ee phi_Sun -
% ??????????????????????????????????????????????????????????????????
% if (0 < phi_Sun) && (phi_Sun < pi)
%     T_Sun_adim = 0;
% end

Eq_old = Eq_i;

% -------------------------------------------------------------------------
% Commented by Federico:
% Note: the output L_i, L_o are always such that:
% L_o > L_i
% Eq_a_o(6) < L_o < Eq_a_o(6) + 2*pi
% -------------------------------------------------------------------------

% Average true anomaly corresponding to L_i and L_o for the eclipse
th_ecl_med = mod( (L_i+L_o)/2 - (kep_in(4)+kep_in(5)), 2*pi);



% -------------------------------------------------------------------------
% Propagation
% -------------------------------------------------------------------------
% L_i is the true longitude of the eclipse entry point and L_p_i is
% the true longitude of the initial point of the perigee thrust
% arc.



% -------------------------------------------------------------------------
% Propagation with J2,J3,J4,J5 drag and third body only, over an entire
% period
% -------------------------------------------------------------------------
% -------------------- J2  ------------------------------------------------
% Propagate for J2 from 0 to 2*pi to avoid eclipses
if ecl_flag
    
    [Eq_end_J2, dt] = AnEquin_all_forward_tang_m([Eq_a_o(6) Eq_a_o(6) + 2*pi], Eq_a_o, ...
        0, 0, 0, ...
        0, 0, 0, ...
        alpha_Sun, beta_Sun, ...
        muadim, geopotential, third_body, R, drag, Earth_flat, ill_flag); % Coasting leg
    
    % For the rest of the propagation do not use J2 up to J5, therefore
    geopotential.J2 = 0;
    geopotential.J3 = 0;
    geopotential.J4 = 0;
    geopotential.J5 = 0;
    
    % Do not use drag and third body
    drag.CD = 0;
    third_body.flag_Moon = 0;
    third_body.flag_Sun = 0;
    
else
    Eq_end_J2 = Eq_a_o;
end






% If the eclipse ENTRY point is before the end of the post-apogee coasting
% arc (that is, before the starting of the perigee thrust arc):
if ~isnan(L_i) && (L_i <= L_p_i) 
    
    % If the eclipse EXIT point L_o is after the end of the perigee coasting arc,
    % propagate the coasting up to the ENTRY point and then propagate the eclipse arc.
    % No thrusting during eclipse. 
    if (L_o > L_p_i) 
        % -----------------------------------------------------------------
        % Commented by Federico:
        %             fprintf('a')
        %             [Eq_a_o(6) L_p_i L_i L_o]
        %             [Eq_a_o(6) L_i L_o]
        % -----------------------------------------------------------------
        
        % Store pericenter thrust arc beginning point
        tmp = L_p_i;
        
        % This will be used as final point of the propagation for the
        % eclipse
        L_p_i  = min([L_o L_p_i + dL_p_i*2]);
        
        % Pericenter thrust arc not in eclipse
        dL_p_i = max([0 (tmp + dL_p_i*2 - L_p_i)/2]);
        
        if L_i> Eq_a_o(6)
            
            % Input to AnEquin_all_forward_tang_m for J2 propagation with
            % solar radiation pressure
            % [Eq_a_o(6) L_i]: vector defining the true longitude values
            % for the propagation
            % Eq_a_o: initial point for propagation
            % epsilon_t   = 0
            % beta_t      = 0
            % epsilon_rth = 0
            % beta_rth    = 0
            % alpha_rth   = 0
            % epsilon_in  = T_Sun_adim/m
            % alpha_in    = alpha_Sun
            % beta_in     = beta_Sun
            % muadim: adimensional gravitational parameter
            % J2
            
            % PROPAGATION--------------------------------------------------
            % From initial point Eq_a_o to initial point of the eclipse
            % (L_i) with J2 and SRP
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Before, now commented:
            %    [Eq_a_o, dt] = AnEquin_all_forward_tang_m([Eq_a_o(6) L_i], Eq_a_o, 0, 0, 0, 0, 0, T_Sun_adim/m, alpha_Sun, beta_Sun, muadim, J2, R, drag, ill_flag); % Coasting leg
            % After:
            % Adimensional acceleration of the low thrust engine
            epsilon_t   = thrust.thrust_t / m;
            epsilon_rth = thrust.thrust_rth / m;
            
            % Propagation
            [Eq_a_o, dt] = AnEquin_all_forward_tang_m([Eq_a_o(6) L_i], Eq_a_o, epsilon_t, thrust.beta_t, epsilon_rth, thrust.alpha_rth, thrust.beta_rth, ...
                                                       T_Sun_adim/m, alpha_Sun, beta_Sun, muadim, geopotential, third_body, R, drag, Earth_flat, ill_flag); % Coasting leg
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update mass:
            m = m - abs(thrust.thrust_t + thrust.thrust_rth) / (thrust.Isp * g0) * dt(end);
%             m = m * exp( -(epsilon_t + epsilon_rth) * dt * m_rate);
            % Update time:
            DT = DT + dt;
        end
        
        % PROPAGATION--------------------------------------------------
        % From initial point Eq_a_o (final point of the previous propagation, 
        % if applicable, to final point of the eclipse (L_p_i) - only J2,
        % no SRP, because we are in eclipse!
        [Eq_p_i,dt] = AnEquin_all_forward_tang_m([Eq_a_o(6) L_p_i], Eq_a_o, 0, 0, 0, 0, 0, 0, 0, 0, muadim, geopotential, third_body, R, drag, Earth_flat, ill_flag); % Eclipse leg
        
        % Update time
        DT=DT+dt;
        
        % Update time in eclipse
        t_ecl=t_ecl+t_ecl_switch*dt;

    % If the eclipse exit point is after the end of the post-apogee coasting arc, propagate the coasting, then the eclipse, and finally another coasting arc.
    else 
        % -----------------------------------------------------------------
        % Commented by Federico:
        %             fprintf('b')
        %             ['Eq_a_o(6)' 'L_p_i' 'L_i' 'L_o']
        %             [Eq_a_o(6) L_p_i L_i L_o]
        % -----------------------------------------------------------------
        
        
        % PROPAGATION------------------------------------------------------
        % Eclipse starting point is after the initial true longitude ->
        % Propagate from initial true longitude to point of eclipse
        % starting point (L_i) with J2 and SRP
        if L_i > Eq_a_o(6)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Before, now commented:
            %             [Eq_a_o,dt] = AnEquin_all_forward_tang_m([Eq_a_o(6) L_i], Eq_a_o, 0, 0, 0, 0, 0, T_Sun_adim/m, alpha_Sun, beta_Sun, muadim, J2, R, drag, ill_flag); % 1st Coasting leg
            % After:
            % Adimensional acceleration of the low thrust engine
            epsilon_t   = thrust.thrust_t / m;
            epsilon_rth = thrust.thrust_rth / m;
            [Eq_a_o,dt] = AnEquin_all_forward_tang_m([Eq_a_o(6) L_i], Eq_a_o, epsilon_t, thrust.beta_t, epsilon_rth, thrust.alpha_rth, thrust.beta_rth, ...
                                                      T_Sun_adim/m, alpha_Sun, beta_Sun, muadim, geopotential, third_body, R, drag, Earth_flat, ill_flag); % 1st Coasting leg
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            DT = DT + dt;
            % Update mass:
%             m = m * exp( -(epsilon_t + epsilon_rth) * dt * m_rate);
            m = m - abs(thrust.thrust_t + thrust.thrust_rth) / (thrust.Isp * g0) * dt(end);
        end
        
        % PROPAGATION------------------------------------------------------
        % Propagate up to final point of eclipse with J2 but no SRP
        % (because we are in the eclipse)! 
        [Eq_a_o,dt] = AnEquin_all_forward_tang_m([Eq_a_o(6) L_o], Eq_a_o, 0, 0, 0, 0, 0, 0, 0, 0, muadim, geopotential, third_body, R, drag, Earth_flat, ill_flag); % Eclipse leg
        DT = DT + dt;
        
        % Se l'eclisse e' stata superata
        if L_i < Eq_a_o(6)
             % Current total eclipse time has already been accounted for, don't count it anymore.
            t_ecl_switch = 0;          
            t_ecl = t_ecl + dt_ecl;
        else
            t_ecl = t_ecl + dt;
        end
        
        if ecl_flag
            % Update Sun angles alpha and beta for new position of the
            % spacecraft
            [alpha_Sun, beta_Sun, Eq_old] = upd_alpha_beta_Sun(alpha_Sun, beta_Sun, Eq_a_o, Eq_old, muadim);
        end
        
        % PROPAGATION------------------------------------------------------
        % Propagate from final point of eclipse condition up to initial
        % point of perigee thrust arc (identified by L_p_i)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Before, now commented:
        %         [Eq_p_i,dt] = AnEquin_all_forward_tang_m([Eq_a_o(6) L_p_i], Eq_a_o, 0, 0, 0, 0, 0, T_Sun_adim/m, alpha_Sun, beta_Sun, muadim, J2, R, drag, ill_flag); % 2nd Coasting leg
        % After:
        % Adimensional acceleration of the low thrust engine
        epsilon_t   = thrust.thrust_t / m;
        epsilon_rth = thrust.thrust_rth / m;
        [Eq_p_i,dt] = AnEquin_all_forward_tang_m([Eq_a_o(6) L_p_i], Eq_a_o, epsilon_t, thrust.beta_t, epsilon_rth, thrust.alpha_rth, thrust.beta_rth, ...
                                                 T_Sun_adim/m, alpha_Sun, beta_Sun, muadim, geopotential, third_body, R, drag, Earth_flat, ill_flag); % 2nd Coasting leg
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DT = DT + dt;
        % Update mass:
%         m = m * exp( -(epsilon_t + epsilon_rth) * dt * m_rate);
        m = m - abs(thrust.thrust_t + thrust.thrust_rth) / (thrust.Isp * g0) * dt(end);
    end
    
else
    % This is the situation in which the eclipse ENTRY point is after the end of the post-apogee coasting
    % arc (is after the beginning of the perigee thrust arc)
    if abs(L_p_i-Eq_a_o(6))>1e-8
        
        % PROPAGATION------------------------------------------------------
        % Propagate from initial point Eq_a_o to initial point of perigee
        % thrust arc (L_p_i) with J2 e SRP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Before, now commented:
        %         [Eq_p_i,dt] = AnEquin_all_forward_tang_m([Eq_a_o(6) L_p_i],Eq_a_o,0,0,0,0,0,T_Sun_adim/m,alpha_Sun,beta_Sun,muadim,J2,R,drag,ill_flag);
        % After
        % Adimensional acceleration of the low thrust engine
        epsilon_t   = thrust.thrust_t / m;
        epsilon_rth = thrust.thrust_rth / m;
        [Eq_p_i,dt] = AnEquin_all_forward_tang_m([Eq_a_o(6) L_p_i],Eq_a_o, epsilon_t, thrust.beta_t, epsilon_rth, thrust.alpha_rth, thrust.beta_rth, ...
                                                T_Sun_adim/m,alpha_Sun,beta_Sun,muadim, geopotential, third_body, R, drag,Earth_flat,  ill_flag);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DT=DT+dt;
        % Update mass:
%         m = m * exp( -(epsilon_t + epsilon_rth) * dt * m_rate);
        m = m - abs(thrust.thrust_t + thrust.thrust_rth) / (thrust.Isp * g0) * dt(end);
    else
        Eq_p_i = Eq_a_o;
    end
end

%% Thrusting at pericenter



% Adimensional acceleration of the low thrust engine
epsilon = T_adim / m;


% ========================== ECLIPSE ======================================
% Check shadow position again and eventually adjust beginning of thrusting arc

% If true longitude variation of perigee thrusting arc is not zero:
if dL_p_i > 0
    
    % ======================== capire =====================================
    % -----------------------------------------------------------------
    % Commented by Federico:
    % [L_i,L_o]=sol_Eclipse_Eq(Eq_p_i,t_0+(ttt+DT)*t/86400,muadim,R);
    % -----------------------------------------------------------------
    if L_o - Eq_p_i(6)<=0
        while L_o-Eq_p_i(6)<=0
            L_i=L_i+2*pi;
            L_o=L_o+2*pi;
        end
    else
        while L_o-Eq_p_i(6)>2*pi
            L_i=L_i-2*pi;
            L_o=L_o-2*pi;
        end
    end
    
    
    % If there is an eclipse
    if ~isnan(L_i)
        
        % ------> L_i --------> Eq_p_i(6) ----------> ????
        % Following condition means: the ENTRY in the eclipse L_i is before the
        % perigee thrusting arc initial point (Eq_p_i(6)) and the EXIT
        % point from the eclipse (L_o) is after the perigee thrusting arc
        % initial point (Eq_p_i(6)) but we don't know yet if the EXIT point
        % is inside the thrusting arc of after the end of the thrusting arc
        if (L_i < Eq_p_i(6)) && (L_o > Eq_p_i(6))
            
            
            % ------> L_i --------> Eq_p_i(6) ---------> L_o -----------> Eq_p_i(6) + dL_p_i * 2 ----------------
            % Case 1p: thrusting starting point is in eclipse, eclipse EXIT point L_o is still within the thrusting arc. 
            % Action: propagate another coasting arc (because we do not thrust while in eclipse) and adjust starting point.
            if L_o < (Eq_p_i(6) + dL_p_i*2)
                % -----------------------------------------------------
                % Commented by Federico:
                %                 fprintf('c')
                %                 [Eq_p_i(6) L_i L_o]
                % -----------------------------------------------------

                % =========== Propagate coasting arc until eclipse EXIT ===
                % Adjust thursting arc duration
                dL_p_i = (Eq_p_i(6) + dL_p_i*2 - L_o) / 2;
                
                % PROPAGATION----------------------------------------------
                % Propagation from Eq_p_i(6) up to L_o (final point of the
                % eclipse) with coast arc (all epsilon put to zero in
                % AnEquin_all_forward_tang_m) and no SRP
                [Eq_p_i, dt] = AnEquin_all_forward_tang_m([Eq_p_i(6) L_o], Eq_p_i, 0, 0, 0, 0, 0, 0, 0, 0, muadim, geopotential, third_body, R, drag, Earth_flat, ill_flag);
                DT = DT + dt;
                t_ecl = t_ecl + t_ecl_switch * dt;
                % =========================================================

                % Check the position of the next eclipse (necessary for
                % long eclipses and long thrusting arcs)
                % [L_i,L_o]=sol_Eclipse_Eq([Eq_p_i(1:5);Eq_p_i(6)+2*dL_p_i],t_0+(ttt+DT)*t/86400,muadim,R);
                while L_o-Eq_p_i(6)+2*dL_p_i<=0
                    L_i=L_i+2*pi;
                    L_o=L_o+2*pi;
                end
                while L_o-Eq_p_i(6)+2*dL_p_i>2*pi
                    L_i=L_i-2*pi;
                    L_o=L_o-2*pi;
                end
                if L_i<=Eq_p_i(6)+2*dL_p_i
%                     fprintf('*')
                   dL_p_i=(L_i-Eq_p_i(6))/2;
                end
            
            % ------> Eq_p_i(6) -------> L_i ---------> L_o ------------> Eq_p_i(6) + dL_p_i * 2 --------
            % Case 2p: starting point of thrust arc is in eclipse, exit point of the eclipse is outside the thrusting arc. 
            % Action: skip the thrusting arc altogether.
            % Non propago con J2?!?!?
            else 
%                 fprintf('d')
                dL_p_i=0;
            end
        
        % ------> L_i --------> Eq_p_i(6) -------> L_o ----------> Eq_p_i(6) + dL_p_i * 2 ------
        % Case 3p: starting point is not eclipse, exit point of the eclipse is within the thrusting arc. 
        % Action: a thrusting and coasting arc (then another thrusting, see later)  
        elseif (L_i < Eq_p_i(6) + dL_p_i*2) && (L_o < Eq_p_i(6) + dL_p_i*2) 
%             fprintf('e')
            dL_p_i=(Eq_p_i(6)+dL_p_i*2-L_o)/2;
            n_FPET=ceil((L_i-Eq_p_i(6))/dL_subarc);
            if ecl_flag
                [alpha_Sun,beta_Sun,Eq_old]=upd_alpha_beta_Sun(alpha_Sun,beta_Sun,Eq_p_i,Eq_old,muadim);
            end
            
            % PROPAGATION--------------------------------------------------
            % Propagation from Eq_p_i(6) up to L_i (initial point of the
            % eclipse) with thrust arc (different file used for the
            % propagation!!!! why?)
%             keyboard
            [Eq_p_i dt] = AnEquin_all_forward_tang_om_change_in_thrust_multi_m2([Eq_p_i(6) L_i], Eq_p_i, k_p*epsilon, -beta_p_i, u_ratio_i, 0, 0, 0,...
                                                                                T_Sun_adim/m,alpha_Sun,beta_Sun,...
                                                                                muadim,geopotential, third_body, R,drag,Earth_flat, ill_flag,n_FPET);
            DT = DT + dt;
            
            % Update consumed mass
%             m = m * exp(-epsilon * dt * m_rate);
            m = m - abs(thrust.thrust_t + thrust.thrust_rth) / (thrust.Isp * g0) * dt(end);
            
            % PROPAGATION--------------------------------------------------
            % Propagation from Eq_p_i(6) up to L_o (final point of the
            % eclipse) with J2
            [Eq_p_i,dt]=AnEquin_all_forward_tang_m([Eq_p_i(6) L_o],Eq_p_i,0,0,0,0,0,0,0,0,muadim,geopotential, third_body, R,drag,Earth_flat, ill_flag);
%             keyboard
            DT=DT+dt;
            t_ecl=t_ecl+t_ecl_switch*dt;
            
        % ----> Eq_p_i(6) -------> L_i ----------> Eq_p_i(6) + dL_p_i * 2------------> L_o ------   
        % Case 4p: starting point is not eclipse, exit point is outside the thrusting arc. 
        % Action: end the thrusting arc prematurely.    
        elseif ((L_i < Eq_p_i(6) + dL_p_i*2) && ( L_o >= Eq_p_i(6) + dL_p_i*2))%||((L_i+2*pi<Eq_p_i(6)+dL_p_i*2)&&(L_o+2*pi>=Eq_p_i(6)+dL_p_i*2)) 
%             fprintf('f')
%             [Eq_p_i(6) L_i L_o]

            % Propagation below
            
            % Duration of thrust arc is redefined - it ends when the
            % eclipse begins
            dL_p_i = (L_i - Eq_p_i(6)) / 2;

        end
    end
end



% =========================================================================
% Thrusting at pericenter
% =========================================================================

% Propagate using thrust if dL_p_i, semiamplitude of the thrust arc, is not
% greater than 0

if dL_p_i > 0
    
    n_FPET=ceil((dL_p_i*2)/dL_subarc);
    
    if ecl_flag
        [alpha_Sun,beta_Sun,Eq_old]=upd_alpha_beta_Sun(alpha_Sun,beta_Sun,Eq_p_i,Eq_old,muadim);
    end
    
    % PROPAGATION----------------------------------------------------------
    % Propagate up to final point of perigee thrust arc
    % In the following propagation epsilon, the acceleration, is multiplied
    % by k_p to obtain perigee raising or not
    [Eq_p_o dt_p] = AnEquin_all_forward_tang_om_change_in_thrust_multi_m2([Eq_p_i(6) Eq_p_i(6) + dL_p_i*2], Eq_p_i, ...
                                                                          k_p*epsilon, -beta_p_i, u_ratio_i, ...
                                                                          epsilon, alfa_p_i, beta_p_i, ...
                                                                            T_Sun_adim/m, alpha_Sun, beta_Sun, ...
                                                                          muadim, geopotential, third_body, R, drag, Earth_flat, ill_flag, n_FPET);                                                              
    
    % ------------------1st Case ------------------------------------------
    % if final eccentricity  is greater than 0.98
    if sqrt( Eq_p_o(2)^2 + Eq_p_o(3)^2 ) >= 0.98

%         save(['Escape_reached_' num2str(round(100*rand))])

        [dL_i,fval,exitflag] = fsolve(@(x)match_Eq_e(x,Eq_p_i,k_p*epsilon,-beta_p_i,u_ratio_i,T_Sun_adim/m,alpha_Sun,beta_Sun,muadim,geopotential, third_body, R,drag,Earth_flat, ill_flag,0.98,n_FPET),0,options);

        if exitflag<1  
            [Eq_f, dt]=AnEquin_all_forward_tang_om_change_in_thrust_multi_m2([Eq_p_i(6) Eq_p_i(6)+dL_i],Eq_p_i,k_p*epsilon,-beta_p_i,u_ratio_i,0,0,0,T_Sun_adim/m,alpha_Sun,beta_Sun,...
                muadim,geopotential, third_body, R,drag, Earth_flat, ill_flag,n_FPET);
            Eq_f=Eq_p_i;
        else
            [Eq_f, dt]=AnEquin_all_forward_tang_om_change_in_thrust_multi_m2([Eq_p_i(6) Eq_p_i(6)+dL_i],Eq_p_i,k_p*epsilon,-beta_p_i,u_ratio_i,0,0,0,T_Sun_adim/m,alpha_Sun,beta_Sun,...
                muadim,geopotential, third_body, R,drag,Earth_flat, ill_flag,n_FPET);
            dt = dt(end);
            % Updated mass
            m = m*exp(-epsilon*dt*m_rate);
            % Update time
            DT = DT + dt;
        end
        
        dL_in=x(6);
        
        % =================================================================
        % Computation of time within the belt
        % =================================================================
        
        % Mean(?) semimajor axis
        a_m = (Eq_i(1) + Eq_f(1)) / 2;
        
        % Mean (?) eccentricity
        e_m = (norm(Eq_i(2:3)) + norm(Eq_f(2:3))) / 2;
        
        if ( a_m*(1+e_m) ) <= r_belt % All below
            t_belt = DT;
        elseif ( a_m*(1-e_m) ) >= r_belt % All above
            t_belt=0;
        else % Generic
            th_belt=acos((a_m*(1-e_m^2)-r_belt)/(e_m*r_belt));
            t_belt=kepEq_t(th_belt,a_m,e_m,muadim,-th_belt,0);
        end

        % Variation of the vector x over the time so far elapsed (DT)
        % Remember that initial conditions are in the variable x
        dx_dt = ([Eq_f(1:5)-x(1:5); dL_in; m-x(7); t_belt; t_ecl]) / DT;

        % Perche'?
        dx_dt(2:3) = NaN;
        

        if nargout>1
            keyboard
            t_f=t+dt;
            x_f=[Eq_f;m;t_belt+x(8);x(9)+t_ecl];
            % -------------------------------------------------------------
            % Commentato da me finche' non capisco cosa fa data_ecl
%             if data_ecl(1,end)<=ttt
%                 data_ecl=data_ecl(:,data_ecl(1,:)<ttt);
%             end
%             data_ecl=[data_ecl [ttt;t_ecl;th_ecl_med]];
            % -------------------------------------------------------------
        end
        return
        
    % If any element of the final equinoctial elements or time is NaN
    elseif any(isnan([Eq_p_o;dt_p]))
        Eq_f=Eq_p_i;
        e=norm(Eq_f(2:3));
        if e~=0
            Eq_f(2:3)=Eq_f(2:3)/norm(Eq_f(2:3));
        else
            Eq_f(2:3)=[1;0];
        end
        dt_p=pi*sqrt(Eq_p_i(1)^3/muadim)/dL_p_i;
        m=m*exp(-epsilon*dt_p*m_rate);
        DT=DT+dt_p;
        
        dL_in=x(6);
        
        % Computation of time within the belt
        a_m=(Eq_i(1)+Eq_f(1))/2;
        e_m=(norm(Eq_i(2:3))+norm(Eq_f(2:3)))/2;
        if (a_m*(1+e_m))<=r_belt % All below
            t_belt=DT;
        elseif (a_m*(1-e_m))>=r_belt % All above
            t_belt=0;
        else % Generic
            th_belt=acos((a_m*(1-e_m^2)-r_belt)/(e_m*r_belt));
            t_belt=kepEq_t(th_belt,a_m,e_m,muadim,-th_belt,0);
        end

        dx_dt=([Eq_f(1:5)-x(1:5); dL_in; m-x(7); t_belt; t_ecl])/DT;
        return
    end
    dt_p=dt_p(end);
    
% If dL_p_i > 0
else
    Eq_p_o=Eq_p_i;
    dt_p=0;
end

% Update mass based on propellant consumption
m = m * exp(-epsilon * dt_p * m_rate);

% DV=DV+dt_p*epsilon;
DT = DT + dt_p;


%% Coasting after pericenter with J2 perturbation

% % Adjust L_a_i such that it is increasing w.r.ttt. Eq_a_o
L_a_i = Eq_i(6) + 2*(pi-dL_a_i);

% -------------------------------------------------------------------------
% Commented by Federico:
% while (L_a_i-Eq_p_o(6))>=2*pi
%     L_a_i=L_a_i-2*pi;
% end
% while (L_a_i-Eq_p_o(6))<0
%     L_a_i=L_a_i+2*pi;
% end
% -------------------------------------------------------------------------
                
% Check shadow position and eventually adjust beginning of thrusting arc
% [L_i,L_o,kep,dt_ecl]=sol_Eclipse_Eq(Eq_p_o,t_0+(ttt+DT)*t/86400,muadim,R);
if L_o-Eq_p_o(6)<=0
    while L_o-Eq_p_o(6)<=0
        L_i=L_i+2*pi;
        L_o=L_o+2*pi;
    end
else
    while L_o-Eq_p_o(6)>2*pi
        L_i=L_i-2*pi;
        L_o=L_o-2*pi;
    end
end



% If the eclipse entry point is before the end of the post-perigee coasting arc
if ~isnan(L_i) && (L_i <= L_a_i) 
%     keyboard
    % ------> L_i --------> L_a_i -------> L_o ----------
    % If the eclipse exit point is after the end of the post-perigee coasting arc, propagate the coasting and then the eclipse arc.
    if(L_o > L_a_i) 
        
        % -----------------------------------------------------------------
        % Commented by Federico:
        %                     fprintf('g')
        %                     ['Eq_p_o(6)' 'L_a_i' 'L_i' 'L_o']
        %                     [Eq_p_o(6) L_a_i L_i L_o]
        % -----------------------------------------------------------------
        tmp = L_a_i;
        L_a_i = min([L_o L_a_i+dL_a_i*2]);
        dL_a_i = (tmp+dL_a_i*2-L_a_i)/2;
                
        if L_i > Eq_p_o(6)
            if ecl_flag
                [alpha_Sun,beta_Sun,Eq_old]=upd_alpha_beta_Sun(alpha_Sun,beta_Sun,Eq_p_o,Eq_old,muadim);
            end
            
            % PROPAGATION--------------------------------------------------
            % Propagate from Eq_p_o to initial point of eclipse, L_i, using
            % J2 and SRP
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Before, now commented:
            %             [Eq_p_o,dt] = AnEquin_all_forward_tang_m([Eq_p_o(6) L_i],Eq_p_o,0,0,0,0,0,T_Sun_adim/m,alpha_Sun,beta_Sun,muadim,J2,R,drag,ill_flag); % Coasting leg
            % After:
            % Adimensional acceleration of the low thrust engine
            epsilon_t   = thrust.thrust_t / m;
            epsilon_rth = thrust.thrust_rth / m;
            [Eq_p_o,dt] = AnEquin_all_forward_tang_m([Eq_p_o(6) L_i], Eq_p_o, epsilon_t, thrust.beta_t, epsilon_rth, thrust.alpha_rth, thrust.beta_rth,...
                                                     T_Sun_adim/m,alpha_Sun,beta_Sun,muadim,geopotential, third_body, R,drag,Earth_flat, ill_flag); % Coasting leg
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            DT = DT+dt;
            % Update mass:
%             m = m * exp( -(epsilon_t + epsilon_rth) * dt * m_rate);
            m = m - abs(thrust.thrust_t + thrust.thrust_rth) / (thrust.Isp * g0) * dt(end);
        end
        
        % PROPAGATION------------------------------------------------------
        % Propagate from Eq_p_p to initial point of apogee thrusting arc,
        % L_a_i, using only J2 because from L_i we are in eclipse#
        [Eq_a_i,dt] = AnEquin_all_forward_tang_m([Eq_p_o(6) L_a_i],Eq_p_o,0,0,0,0,0,0,0,0,muadim,geopotential, third_body, R,drag,Earth_flat, ill_flag); % Eclipse leg

        DT = DT + dt;
        
        t_ecl = t_ecl + t_ecl_switch*dt;

    % ------> L_i --------> L_o -------> L_a_i ----------
    % If the eclipse exit point is before the begin of the apogee thrust arc, propagate the coasting, then the eclipse, and finally another coasting arc.
    else 
        % -----------------------------------------------------------------
        % Commented by Federico:
        %                     fprintf('h')
        %             [Eq_p_o(6) L_i L_o]
        % -----------------------------------------------------------------
        
        if L_i > Eq_p_o(6)
            if ecl_flag
                [alpha_Sun,beta_Sun,Eq_old]=upd_alpha_beta_Sun(alpha_Sun,beta_Sun,Eq_p_o,Eq_old,muadim);
            end
            
            % PROPAGATION--------------------------------------------------
            % Propagate from Eq_p_o to initial point of the eclipse L_i
            % with J2 and SRP
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Before, now commented:
            %             [Eq_p_o,dt]=AnEquin_all_forward_tang_m([Eq_p_o(6) L_i],Eq_p_o,0,0,0,0,0,T_Sun_adim/m,alpha_Sun,beta_Sun,muadim,J2,R,drag,ill_flag); % 1st Coasting leg
            % After:
            % Adimensional acceleration of the low thrust engine
%             keyboard
            epsilon_t   = thrust.thrust_t / m;
            epsilon_rth = thrust.thrust_rth / m;
            [Eq_p_o,dt]=AnEquin_all_forward_tang_m([Eq_p_o(6) L_i],Eq_p_o, epsilon_t, thrust.beta_t, epsilon_rth, thrust.alpha_rth, thrust.beta_rth,...
                                                    T_Sun_adim/m,alpha_Sun,beta_Sun,muadim,geopotential, third_body, R,drag,Earth_flat, ill_flag); % 1st Coasting leg
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            DT=DT+dt;
            % Update mass:
%             m = m * exp( -(epsilon_t + epsilon_rth) * dt * m_rate);
            m = m - abs(thrust.thrust_t + thrust.thrust_rth) / (thrust.Isp * g0) * dt(end);
        end
        
        % PROPAGATION------------------------------------------------------
        % Propagate from Eq_p_o to final point of the eclipse L_o with J2
        % and no SRP because we are in the eclipse
        [Eq_p_o,dt]=AnEquin_all_forward_tang_m([Eq_p_o(6) L_o],Eq_p_o,0,0,0,0,0,0,0,0,muadim,geopotential, third_body, R,drag,Earth_flat, ill_flag); % Eclipse leg

        DT=DT+dt;
        t_ecl=t_ecl+t_ecl_switch*dt;
%         keyboard
        if ecl_flag
            [alpha_Sun,beta_Sun,Eq_old]=upd_alpha_beta_Sun(alpha_Sun,beta_Sun,Eq_p_o,Eq_old,muadim);
        end
        
        % PROPAGATION------------------------------------------------------
        % Propagate from Eq_p_o to initial point L_a_i of the apogee
        % thrusting arc using J2 and SRP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Before, now commented:
        %         [Eq_a_i,dt]=AnEquin_all_forward_tang_m([Eq_p_o(6) L_a_i],Eq_p_o,0,0,0,0,0,T_Sun_adim/m,alpha_Sun,beta_Sun,muadim,J2,R,drag,ill_flag); % 2nd Coasting leg
        % After:
        % Adimensional acceleration of the low thrust engine
        epsilon_t   = thrust.thrust_t / m;
        epsilon_rth = thrust.thrust_rth / m;
        [Eq_a_i,dt]=AnEquin_all_forward_tang_m([Eq_p_o(6) L_a_i],Eq_p_o, epsilon_t, thrust.beta_t, epsilon_rth, thrust.alpha_rth, thrust.beta_rth,...
                                               T_Sun_adim/m,alpha_Sun,beta_Sun,muadim,geopotential, third_body, R,drag,Earth_flat, ill_flag); % 2nd Coasting leg
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DT=DT+dt;
        % Update mass:
%         m = m * exp( -(epsilon_t + epsilon_rth) * dt * m_rate);
        m = m - abs(thrust.thrust_t + thrust.thrust_rth) / (thrust.Isp * g0) * dt(end);
    end
  
% If the eclipse entry point is AFTER the end of the post-perigee coasting arc
else
    
    if abs(L_a_i-Eq_p_o(6))>1e-8
        if ecl_flag
            [alpha_Sun,beta_Sun,Eq_old]=upd_alpha_beta_Sun(alpha_Sun,beta_Sun,Eq_p_o,Eq_old,muadim);
        end
        
        % PROPAGATION------------------------------------------------------
        % Propagate from Eq_p_o to initial point L_a_i of the apogee
        % thrusting arc using J2 and SRP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Before, now commented:
        %         [Eq_a_i,dt]=AnEquin_all_forward_tang_m([Eq_p_o(6) L_a_i],Eq_p_o,0,0,0,0,0,T_Sun_adim/m,alpha_Sun,beta_Sun,muadim,J2,R,drag,ill_flag);
        % After:
        % Adimensional acceleration of the low thrust engine
        epsilon_t   = thrust.thrust_t / m;
        epsilon_rth = thrust.thrust_rth / m;
        [Eq_a_i,dt]=AnEquin_all_forward_tang_m([Eq_p_o(6) L_a_i],Eq_p_o, epsilon_t, thrust.beta_t, epsilon_rth, thrust.alpha_rth, thrust.beta_rth,...
                                               T_Sun_adim/m,alpha_Sun,beta_Sun,muadim,geopotential, third_body, R,drag,Earth_flat, ill_flag);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         keyboard
        DT=DT+dt;
        % Update mass:
%         m = m * exp( -(epsilon_t + epsilon_rth) * dt * m_rate);
        m = m - abs(thrust.thrust_t + thrust.thrust_rth) / (thrust.Isp * g0) * dt(end);
    else
        Eq_a_i=Eq_p_o;
    end
end




%% Thrusting at apocenter

epsilon = T_adim/m;


% Eventually adjust thrusting arc

% If thrusting arc at apogee is not zero
if dL_a_i>0
    
    if L_o-Eq_a_i(6)<=0
        while L_o-Eq_a_i(6)<=0
            L_i=L_i+2*pi;
            L_o=L_o+2*pi;
        end
    else
        while L_o-Eq_a_i(6)>2*pi
            L_i=L_i-2*pi;
            L_o=L_o-2*pi;
        end
    end
    
    % If there is an eclipse
    if ~isnan(L_i)
        
        % --------- L_i -----> Eq_a_i(6) ---------
        % Initial point of apogee thrusting arc is after initial eclipse
        % point but we don't know yet if the eclipse ends before or after
        % the end of the thrusting arc
        if (L_i < Eq_a_i(6)) && (L_o > Eq_a_i(6))
            
            
            % --------- L_i -----> Eq_a_i(6) ---------> L_o ------------> Eq_a_i(6) + 2 * dL_a_i --------------
            % Case 1a: starting point of the thrusting arc is in eclipse, exit point of the eclipse is within the thrusting arc.
            % Action: propagate another coasting arc (until L_o) and adjust starting point.
            if (L_o < Eq_a_i(6) + dL_a_i*2) 
                
                % ---------------------------------------------------------
                % Commented by Federico:
                %                 fprintf('i')
                %                 [Eq_a_i(6) L_i L_o]
                % ---------------------------------------------------------
                % Adjust length of the thrusting arc
                dL_a_i = (Eq_a_i(6) + dL_a_i*2 - L_o)/2;
                
                % PROPAGATION------------------------------------------------------
                % Propagate from Eq_a_i to final point of the eclipse L_o
                % using only J2 (no SRP, eclipse!)
                [Eq_a_i,dt]=AnEquin_all_forward_tang_m([Eq_a_i(6) L_o],Eq_a_i,0,0,0,0,0,0,0,0,muadim,geopotential,third_body,R,drag,Earth_flat, ill_flag);
                keyboard
                DT=DT+dt;
                
                t_ecl=t_ecl+t_ecl_switch*dt;
%                 keyboard
                % Check the position of the next eclipse (necessary for
                % long eclipses and long thrusting arcs)
                % [L_i,L_o]=sol_Eclipse_Eq([Eq_a_i(1:5);Eq_a_i(6)+2*dL_a_i],t_0+(ttt+DT)*t/86400,muadim,R);
                
                while L_o-Eq_a_i(6)+2*dL_p_i<=0
                    L_i=L_i+2*pi;
                    L_o=L_o+2*pi;
                end
                while L_o-Eq_a_i(6)+2*dL_p_i>2*pi
                    L_i=L_i-2*pi;
                    L_o=L_o-2*pi;
                end
                if L_i<=Eq_a_i(6)+2*dL_a_i
%                     fprintf('*')
                    L_adj=Eq_a_i(6)+2*dL_a_i;
                    dL_a_i=(L_i-Eq_a_i(6))/2;
                    adj_flag=1;
%                     keyboard
                end
            
            % --------- L_i -----> Eq_a_i(6) ---------> Eq_a_i(6) + 2 * dL_a_i -------------> L_o ----------
            % Case 2a: starting point is in eclipse, exit point from eclipse is outside the thrusting arc. 
            % Action: skip the thrusting arc altogether and propagate a coasting arc instead.
            else 
%                 fprintf('j')

                % PROPAGATION----------------------------------------------
                % Propagate from Eq_a_i to final point of the thrusting arc
                % with no thrust! but only with J2 (no SRP because we are
                % in eclipse)
                [Eq_a_i,dt]=AnEquin_all_forward_tang_m([Eq_a_i(6) Eq_a_i(6)+2*dL_a_i],Eq_a_i,0,0,0,0,0,0,0,0,muadim,geopotential,third_body,R,drag,Earth_flat, ill_flag);
                keyboard
                DT=DT+dt;
                
                t_ecl=t_ecl+t_ecl_switch*dt;
                dL_a_i=0;
            end
           
        % --------- Eq_a_i(6) -----> L_i ---------> L_o-------------> Eq_a_i(6) + 2 * dL_a_i  ----------
        % Case 3a: starting point is not eclipse, exit point from eclipse is within the thrusting arc. 
        % Action: a thrusting and coasting arc.
        elseif (L_i < Eq_a_i(6) + dL_a_i*2) && (L_o < Eq_a_i(6) + dL_a_i*2) 
            
            % -------------------------------------------------------------
            % Commented by Federico:
            %             fprintf('k')
            %             [Eq_a_i(6) Eq_a_i(6)+2*dL_a_i L_i L_o]
            % -------------------------------------------------------------
            
            % Redefine dL_a_i for the thrust arc to follow
            dL_a_i = (Eq_a_i(6)+dL_a_i*2-L_o)/2;
            n_FPET=ceil((L_i-Eq_a_i(6))/dL_subarc);
            if ecl_flag
                [alpha_Sun,beta_Sun,Eq_old]=upd_alpha_beta_Sun(alpha_Sun,beta_Sun,Eq_a_i,Eq_old,muadim);
            end
            
            % PROPAGATION--------------------------------------------------
            % Propagate from Eq_a_i to initial point of eclipse, L_i, with
            % J2 and SRP and Thrust!
            [Eq_a_i dt]=AnEquin_all_forward_tang_om_change_in_thrust_multi_m2([Eq_a_i(6) L_i],Eq_a_i,k_a*epsilon,beta_a_i,u_ratio_i,0,0,0,T_Sun_adim/m,alpha_Sun,beta_Sun,...
                                                                               muadim,geopotential,third_body,R,drag,Earth_flat, ill_flag,n_FPET);
%             keyboard
            DT=DT+dt;
            
            % Update mass according to mass consumption
%             m=m*exp(-epsilon*dt*m_rate);
            m = m - abs(thrust.thrust_t + thrust.thrust_rth) / (thrust.Isp * g0) * dt(end);
            
            % PROPAGATION--------------------------------------------------
            % Propagate from Eq_a_i to final point of eclipse L_o using
            % only J2 (we are in the eclipse!)
            [Eq_a_i,dt]=AnEquin_all_forward_tang_m([Eq_a_i(6) L_o],Eq_a_i,0,0,0,0,0,0,0,0,muadim,geopotential,third_body,R,drag,Earth_flat, ill_flag);
            keyboard
            DT=DT+dt;
            t_ecl=t_ecl+t_ecl_switch*dt;
            
         % --------- Eq_a_i(6) -----> L_i ---------> Eq_a_i(6) + 2 * dL_a_i -------------> L_o ----------
         % Case 4a: starting point is not eclipse, exit point from eclipse is outside the thrusting arc. 
         % Action: end the thrusting arc prematurely and propagate a coasting arc afterwards.
        elseif ((L_i < Eq_a_i(6) + dL_a_i*2) && (L_o >= Eq_a_i(6) + dL_a_i*2))%||((L_i+2*pi<Eq_a_i(6)+dL_a_i*2)&&(L_o+2*pi>=Eq_a_i(6)+dL_a_i*2))
            
            % -------------------------------------------------------------
            % Commented by Federico:
            %             fprintf('l')
            %             [Eq_a_i(6) Eq_a_i(6)+2*dL_a_i L_i L_o]
            % -------------------------------------------------------------
            % Redefine semiamplitude of thrust arc at apogee
            L_adj = Eq_a_i(6) + 2*dL_a_i;
            dL_a_i = (L_i-Eq_a_i(6))/2;
            adj_flag=1;
        end
    end
end

% =========================================================================
% Thrusting at apocenter
% =========================================================================
if dL_a_i>0
    
    n_FPET=ceil((dL_a_i*2)/dL_subarc);
    if ecl_flag
        [alpha_Sun,beta_Sun,Eq_old]=upd_alpha_beta_Sun(alpha_Sun,beta_Sun,Eq_a_i,Eq_old,muadim);
    end
    
    % PROPAGATION ---------------------------------------------------------
    % Propagation over the thrust arc of apogee (with amplitude as defined
    % above) using J2, SRP and thrust
    [Eq_f, dt_f]=AnEquin_all_forward_tang_om_change_in_thrust_multi_m2([Eq_a_i(6) Eq_a_i(6)+dL_a_i*2],Eq_a_i,...
                                                                        k_a*epsilon,beta_a_i,u_ratio_i,...
                                                                        epsilon,alfa_a_i,beta_a_i,...
                                                                        T_Sun_adim/m,alpha_Sun,beta_Sun,...
                                                                        muadim,geopotential,third_body,R,drag,Earth_flat, ill_flag,n_FPET);
%     keyboard
    
    % If eccentricity of final orbit is greater than 0.98
    if ((Eq_f(2)^2+Eq_f(3)^2)>=0.98)
        
        [dL_i,fval,exitflag]=fsolve(@(x)match_Eq_e(x,Eq_a_i,k_a*epsilon,beta_a_i,u_ratio_i,T_Sun_adim/m,alpha_Sun,beta_Sun,muadim,geopotential,third_body,R,drag,Earth_flat, ill_flag,0.98,n_FPET),0,options);
        
        if exitflag<1
            fprintf('<')
             Eq_f=Eq_a_i;
        else
            [Eq_f, dt]=AnEquin_all_forward_tang_om_change_in_thrust_multi_m2([Eq_a_i(6) Eq_a_i(6)+dL_i],Eq_a_i,k_a*epsilon,beta_a_i,u_ratio_i,0,0,0,T_Sun_adim/m,alpha_Sun,beta_Sun,...
                                                                             muadim,geopotential,third_body,R,drag,Earth_flat, ill_flag,n_FPET);
            dt=dt(end);
            m=m*exp(-epsilon*dt*m_rate);
%             DV=DV+dt*epsilon;
            DT=DT+dt;
        end
        
        dL_in=x(6);
        
        % =================================================================
        % Computation of time within the belt
        % =================================================================
        a_m=(Eq_i(1)+Eq_f(1))/2;
        e_m=(norm(Eq_i(2:3))+norm(Eq_f(2:3)))/2;
        if (a_m*(1+e_m))<=r_belt % All below
            t_belt=DT;
        elseif (a_m*(1-e_m))>=r_belt % All above
            t_belt=0;
        else % Generic
            th_belt=acos((a_m*(1-e_m^2)-r_belt)/(e_m*r_belt));
            t_belt=kepEq_t(th_belt,a_m,e_m,muadim,-th_belt,0);
        end

        dx_dt=([Eq_f(1:5)-x(1:5);dL_in;m-x(7);t_belt;t_ecl])/DT;
        dx_dt(2:3)=NaN;
        
        if nargout>1
            t_f=t+dt;
            x_f=[Eq_f;m;t_belt+x(8);x(9)+t_ecl];
            % -------------------------------------------------------------
            % Commentato da me finche' non capisco cosa fa data_ecl
%             if data_ecl(1,end)<=ttt
%                 data_ecl=data_ecl(:,data_ecl(1,:)<ttt);
%             end
%             data_ecl=[data_ecl [ttt;t_ecl;th_ecl_med]];
            % -------------------------------------------------------------
            
        end
        return
        
    elseif any(isnan([Eq_f;dt_f]))
        Eq_f=Eq_a_i;
        e=norm(Eq_f(2:3));
        if e~=0
            Eq_f(2:3)=Eq_f(2:3)/norm(Eq_f(2:3));
        else
            Eq_f(2:3)=[1;0];
        end
        dt_f=pi*sqrt(Eq_f(1)^3/muadim)/dL_p_i;
        m=m*exp(-epsilon*dt_f*m_rate);
%         DV=DV+dt_f*epsilon;
        DT=DT+dt_f;
        
        dL_in=x(6);
        
        % Computation of time within the belt
        a_m=(Eq_i(1)+Eq_f(1))/2;
        e_m=(norm(Eq_i(2:3))+norm(Eq_f(2:3)))/2;
        if (a_m*(1+e_m))<=r_belt % All below
            t_belt=DT;
        elseif (a_m*(1-e_m))>=r_belt % All above
            t_belt=0;
        else % Generic
            th_belt=acos((a_m*(1-e_m^2)-r_belt)/(e_m*r_belt));
            t_belt=kepEq_t(th_belt,a_m,e_m,muadim,-th_belt,0);
        end

        dx_dt=([Eq_f(1:5)-x(1:5);dL_in;m-x(7);t_belt;t_ecl])/DT;
        return
    end
    
    dt_f=dt_f(end);
    
else
    
    Eq_f=Eq_a_i;
    dt_f=0;
end


% Update mass
m = m*exp(-epsilon*dt_f*m_rate);
% DV=DV+dt_f*epsilon;
% Update elapsed time
DT=DT+dt_f;

if adj_flag
%     fprintf('i')
%     [Eq_f(6) L_adj]
    [Eq_f,dt]=AnEquin_all_forward_tang_m([Eq_f(6) L_adj],Eq_f,0,0,0,0,0,0,0,0,muadim,geopotential,third_body,R,drag,Earth_flat, ill_flag);
%     keyboard
    DT=DT+dt;
    t_ecl=t_ecl+t_ecl_switch*dt;
    % -------------------------------------------------------------
    % Commented by Federico:
    %     kep_f=eq2kep(Eq_f);
    %     dL_in2=mod((kep_f(4)+kep_f(5))-(kep_in(4)+kep_in(5)),2*pi);
    %     if dL_in2>3*pi/2
    %         dL_in2=dL_in2-2*pi;
    %     elseif dL_in2>pi/2
    %         dL_in2=dL_in2-pi;
    %     end
    %     dL_in2
    % -------------------------------------------------------------
end

%% Computation of displacement of apses line 
% To be put in the 6th row of the dx_dt variable

% Keplerian elements of everything computed
kep_f = eq2kep(Eq_f);

kep_f_J2 = eq2kep(Eq_end_J2);

% Displacement of apses line
dL_in = mod( (kep_f(4) + kep_f(5)) + (kep_f_J2(4) + kep_f_J2(5)) - 2*(kep_in(4)+kep_in(5) ), 2*pi);

% If initial eccentricity was approximately zero (no omega defined?)
if sqrt(x(2)^2+x(3)^2)<1e-5
    dL_in = 0;
% I think that he brings dL_in in the interval -pi/2 pi/2 here
elseif dL_in>3*pi/2
    dL_in=dL_in-2*pi;
elseif dL_in>pi/2
    dL_in=dL_in-pi;
end

%% Computation of time within the belt
% 8th row of dx_dt

% Average semimajor axis
a_m = (Eq_i(1)+Eq_f(1)) / 2;
% Average eccentricity
e_m = (norm(Eq_i(2:3))+norm(Eq_f(2:3)))/2;

% All below
% If average apogee height is lower than radiation belt, time spent under
% the radiation belt (t_belt) correspond to the total time DT
if (a_m*(1+e_m))<=r_belt 
    t_belt=DT;
% All above
% If average perigee height is greater than radiation belt radius, the time
% spent below the radiation belt radius is zero
elseif (a_m*(1-e_m))>=r_belt 
    t_belt=0;
% Generic
% If none of the above applies -cosa fa qui?
else 
%     keyboard
%     th_belt=acos((a_m*(1-e_m^2)-r_belt)/(e_m*r_belt));
%     t_belt=kepEq_t(th_belt,a_m,e_m,muadim,-th_belt,0);
t_belt = 0;
end


%%
% This is the only place where dx_dt is computed unless particular
% situations verifies at the perigee or apogee thrusting (See above).
% This is the averaged vector, as in Equation 39 of the paper "Extended
% Analytical Formulas for the Perturbed Keplerian Motion Under a Constant
% Control Acceleration"
% Rows 1-5: variation of the first 5 equinoctial elements over the time DT
%           Eq_f are the final equinoctial elements, Eq_i are the initial
%           equinoctial elements
% Row 6: line of apsides displacement
% Row 7: variation of mass
% Row 8: time within the belt
% Row 9: time of eclipse


% 1. Federico:
% He divide elements variations by DT, which is the sum of the times as
% computed by the analytical propagator

% dx_dt_old=([ Eq_f(1:5,end) - Eq_i(1:5); ...
%             dL_in; ...
%             m - x(7); ...
%             t_belt; ...
%             t_ecl])/DT;
        


% Semimajor axis
a = Eq_i(1);
% Maan motion
n = sqrt(muadim/a^3);
% Eccentricity
e = sqrt(Eq_i(2)^2 + Eq_i(3)^2);
% Parameter
p = a * (1-e^2);
% Inclination
i = Eq_i(3);

% 
n_bar = n * (1 + 1.5 * J2 * (R/p)^2 * sqrt(1-e^2)*(1-1.5*sin(i)*sin(i)));


dx_dt=([ Eq_end_J2(1:5,end) - Eq_a_o(1:5) + Eq_f(1:5,end) - Eq_i(1:5); ...
         dL_in; ...
         m - x(7); ...
         t_belt; ...
         t_ecl])/(2*pi);
% keyboard
  
B = sqrt(1-e^2);
Phi = 1 + Eq_i(3);
dt_dL = sqrt(a^3/muadim) * B^3 / Phi^2;


dx_dt = dx_dt * n_bar;



% -------------------------------------------------------------------------
% Commented by Federico:
% if ttt*t/86400>3696
% % %     figure(99)
% % %     for i=1:7
% % %         subplot(3,3,i)
% % %         plot(ttt*t/86400,x(i))
% % %     end
% % %     subplot(3,3,8)
% % %     plot(ttt*t/86400,th_ecl_med*180/pi)
% % %     subplot(3,3,9)
% % %     plot(ttt*t/86400,(Eq_f(6)-Eq_i(6))*180/pi)
% % %     drawnow
%     figure(100)
% % %     for i=1:9
% % %         subplot(3,3,i)
% % %         plot(ttt*t/86400,dx_dt(i)*DT)
% % %     end
%     plot(ttt*t/86400,dL_in*180/pi)
% %     plot(ttt*t/86400,sqrt(x(2)^2+x(3)^2))
% % %     plot(ttt*t/86400,t_ecl*t/3600)
% % % plot(ttt*t/86400,DT*t/3600)
%     drawnow
% % 
% end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Commentato da me finche' non capisco cosa fa data_ecl
% if data_ecl(1,end)<=ttt
%     data_ecl=data_ecl(:,data_ecl(1,:)<ttt);
% end
% data_ecl=[data_ecl [ttt;t_ecl;th_ecl_med]];
% -------------------------------------------------------------------------

if nargout>1
    t_f=t+DT;
    x_f=[Eq_f;m;t_belt+x(8);x(9)+t_ecl];
end

% -------------------------------------------------------------------------
% Commented by Federico:
% figure(99)
% subplot(1,2,1)
% plot(ttt*t/86400,dx_dt(4),'.')
% subplot(1,2,2)
% plot(ttt*t/86400,dx_dt(5),'.')
% dx_dt
% Eq_i
% Eq_f
% keyboard
% if ttt*t/86400>150
%     keyboard
% end
% [((kep_f(4)+kep_f(5))-(kep_in(4)+kep_in(5)))*180/pi t_ecl]
% fprintf('\n----------\n')
% % -------------------------------------------------------------------------
% keyboard
return







function [alpha_Sun, beta_Sun, Eq_old] = upd_alpha_beta_Sun(alpha_Sun, beta_Sun, Eq, Eq_old, mu)

    % Cartesian elements current equinoctial elements 
    Car = eq2cart(Eq,mu);
    
    % Cartesian elements old equinoctial elements
    Car_old = eq2cart(Eq_old,mu);
    
    % Adimensional Sun position vector in rth reference frame of old Eq_old
    % state vector
    r_tnh_old = [cos(beta_Sun) * cos(alpha_Sun); ...
                 cos(beta_Sun) * sin(alpha_Sun); ...
                 sin(beta_Sun)];
          
    % Adimensional Sun position vector in inertial reference frame         
    % Transformation from rth to inertial reference frame for
    % adimensional Sun position vector
    r_car = rth_carT(r_tnh_old, Car_old);
    
    % Adimensional Sun position vector in rth reference frame of Car state
    % vector
    r_rth = car_rthT(r_car,Car);
    
    % Sun distance from current state position
    rr = norm(r_rth);
    
    % Sun elevation angle with respect to new position vector
    beta_Sun = asin(r_rth(3)/rr);
    
    % Sun azimuth angle wrt new position vector
    alpha_Sun = atan2(r_rth(2)/rr/cos(beta_Sun),r_rth(1)/rr/cos(beta_Sun));
    
    % Update the old equinoctial elements
    Eq_old = Eq;

return






function err=match_Eq_dt(x,Eq_0,epsilon,beta,u_ratio_i,epsilon_Sun,alpha_Sun,beta_Sun,muadim,J2,R,ill_flag,dt_nom,n_FPET)

    [Eq_f, dt]=AnEquin_all_forward_tang_om_change_in_thrust_multi_m([Eq_0(6) Eq_0(6)+x],Eq_0,epsilon,beta,u_ratio_i,0,0,0,epsilon_Sun,alpha_Sun,beta_Sun,muadim,J2,R,ill_flag,n_FPET);
    err=dt-dt_nom;
    
return




function err=match_Eq_e(x,Eq_0,epsilon,beta,u_ratio_i,epsilon_Sun,alpha_Sun,beta_Sun,muadim,geopotential,third_body,R,drag,ill_flag,e_target,n_FPET)

    [Eq_f, dt]=AnEquin_all_forward_tang_om_change_in_thrust_multi_m2([Eq_0(6) Eq_0(6)+x],Eq_0,epsilon,beta,u_ratio_i,0,0,0,epsilon_Sun,alpha_Sun,beta_Sun,muadim,geopotential,third_body, R,drag,ill_flag,n_FPET);
    err=sqrt(Eq_f(2)^2+Eq_f(3)^2)-e_target;
    
return

