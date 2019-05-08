function dEq_dt = equations_motion3(t, x, input, thrust, control, drag_NUM, g0)
% function dEq_dt = equations_motion3(x, input, thrust, drag_NUM, g0)
% Equation of motion: derivative of the equinoctial elements wrt to time 
% Reference: Zuiani et al "Extended Analytical Formulas for..." 
%            Equations [2]

%% Input 


AU =astro_constants(2);
% J2 []
% keyboard
J2 = input.J2;

% Distance unit
DU = input.DU;

% Time unit
TU = input.TU;

% Earth radius []
R = input.R;

% Gravitational acceleration [DU^2/TU^3]
mu = input.muadim;

% Initial date [MJD2000]
t0 = input.t0;

if thrust.epsilon_inertial
    P_Sun = input.SRP.P_Sun;
    Cr = input.SRP.Cr;
    A_m_SRP = input.SRP.A_m;
end



% % Tangential acceleration
% % epsilon_t = input.epsilon_t;

% % Rth acceleration
% epsilon_rth = input.epsilon_rth;
% alpha_rth   = input.alpha_rth;
% beta_rth    = input.beta_rth;
% 
% % Inertial acceleration
% epsilon_In = input.epsilon_In;
% beta_In    = input.beta_In;
% gamma0_In  = input.gamma0_In; 



%% Initialization

% Equinoctial elements
a  = x(1);
P1 = x(2);
P2 = x(3);
Q1 = x(4);
Q2 = x(5);
l  = x(6);

% Mass
m = x(7);

% Eccentricity
e = sqrt(P1^2 + P2^2);

% 
B = sqrt(1 - e^2);

% Compute L from l
% Battin, Eq (10.53)
b = a * sqrt(1 - P1^2 - P2^2);
n = sqrt(mu / a^3);
h = n * a * b;

% Solve Kepler equation in the eccentirc longitude K = Omega + omega + E
K = kepler_K(e, l, P1, P2);

% Radius
r = a * (1 - P1 * sin(K) - P2 * cos(K));

% Compute true longitude L = Omega + omega + theta (Battin, page. 493)
sin_L = a/r * ( (1 - a/(a+b) * P2^2) * sin(K) + a/(a+b) * P1 * P2 * cos(K) - P1 );
cos_L = a/r * ( (1 - a/(a+b) * P1^2) * cos(K) + a/(a+b) * P1 * P2 * sin(K) - P2 );
L = atan2(sin_L,cos_L);


Phi = 1 + P1 * sin(L) + P2 * cos(L);
D = sqrt(1 + P1^2 + P2^2 + 2 * (P1 * sin(L) + P2 * cos(L)));
G = 1 + Q1^2 + Q2^2;

% Parameter [DU]
p = a * (1-e^2);

% Radius [DU]
r = p / Phi;

% True anomaly?
% RAAN:
RAAN_omega = atan2(P1, P2);
theta = mod(L - RAAN_omega, 2*pi);


k_a = control.k_a;
k_p = control.k_p;


%% Eclipse and Sun angles

Eq = [a P1 P2 Q1 Q2 L];

% kep = eq2kep(Eq);
% kep(1) * input.DU  - input.DU
% pause

% [L_i, L_o] = sol_Eclipse_Eq(Eq, t0 + t * TU / 86400, mu, R);
% 
% 
% 
% % Evaluate if spacecraft is in eclipse
% if L >= L_i && L <= L_o
%     eclipse = 1;
% else
%     eclipse = 0;
% end
eclipse = 0;

% 
% % Here I follow what Federico did in sol_Eclipse_Eq to compute the
% % position of the Sun:
% 
% % Earth position vector in heliocentric reference frame [km]
% r_S_E = EphSS(3,t0 + t*TU);
% 
% % Ecliptic inclination [rad]
% incl_ecl = astro_constants(8);
% 
% % Rotate r_S_R: from ecliptic to equatorial reference place. Therefore
% % this is the vector Sun-Earth in the Earth centered reference frame.
% r_S_E = [r_S_E(1); cos(incl_ecl)*r_S_E(2); sin(incl_ecl)*r_S_E(2)];
% 
% % Earth-Sun position vector [km]
% r_E_S = - r_S_E;
% 
% % ---------------------------------------------------------------------
% %     % Otherwise, compute directly sun position (this function is not part
% %     % of the AstroToolBox and has not been officialy validated)
% %     % Transformation from MJD2000 to JD
% %     JD = mjd20002jd(t0 + t*TU);
% %     r_Sun = sun_position(JD)';
% %     % Sun position in Earth reference frame [DU]
% %     r_Sun = r_Sun / DU;
% % ---------------------------------------------------------------------
% 
% % Earth-Sun position vector [DU]
% r_E_S = r_E_S / DU;
% 
% % Spacecraft cartesian coordinates [DU]
% kep_SC = eq2kep(x);
% state_SC = kep2cart(kep_SC,mu);
% 
% % Spacecraft position vector [DU] in Earth reference frame
% r_E_SC = state_SC(1:3);
% 
% % Cylindrical shadow model: Ref: Mengali and Quarta, Fondamenti di
% % Meccanica del Volo Spaziale
% 
% tau = dot(r_E_SC, (r_E_SC' - r_E_S)) /dot((r_E_S - r_E_SC'),(r_E_S - r_E_SC'));
% 
% if tau < 0 || tau > 1
%     % No eclipse
% elseif tau>=0 && tau<=1
%     % Eclipse is possbile
%     c = sqrt(norm(r_E_SC)^2 * (1-tau) + tau * dot(r_E_SC,r_E_S));
%     if c<R
%         eclipse = 1;
%     end
% end
% 
% 
% % Sun angles: from file sol_Eclipse_Eq.m of Federico Zuiani 
% 
% % Sun-Earth vector from Earth reference frame to radial-transversal-h reference frame
% r_S_E_rth = car_rthT(r_S_E, kep2cart(kep_SC,mu));
% 
% norm_r_S_E_rth = norm(r_S_E_rth);
% 
% % Sun elevation angle [rad]
% beta_Sun = asin(r_S_E_rth(3)/norm_r_S_E_rth);
% 
% % Sun azimuth angle [rad]
% alpha_Sun = atan2(r_S_E_rth(2)/norm_r_S_E_rth/cos(beta_Sun),r_S_E_rth(1)/norm_r_S_E_rth/cos(beta_Sun));
%     


%% Accelerations




% =========================================================================
% Accelerations DRAG
% =========================================================================
if drag_NUM.CD 
    % Compute density here!

    % Eccentricity []
    e = sqrt(P1^2 + P2^2);

    % Parameter [DU]
    p = a * (1-e^2);

    % Radius [DU]
    r = p / Phi;

    % Height [DU]
    height = r - R;

    % Height [km]
    height = height * DU;

    % Atmospheric density [kg/m^3]
    rho = exponential_atm_model(height);

    % Atmospheric density [kg/DU^3]
    drag_NUM.rho = (rho/1e-9) * DU^3;

    % Drag accekleration
    aR_Drag = - 0.5 * drag_NUM.rho * drag_NUM.CD * drag_NUM.A_m * ( mu / a ) * (2 * Phi / B^2 - 1) * (P2 * sin(L) - P1 * cos(L)) / D;
    aT_Drag = - 0.5 * drag_NUM.rho * drag_NUM.CD * drag_NUM.A_m * ( mu / a ) * (2 * Phi / B^2 - 1) * (1 + P1 * sin(L) + P2 * cos(L)) / D;
    aN_Drag = 0;
else
    aR_Drag = 0;
    aT_Drag = 0;
    aN_Drag = 0;
end

% =========================================================================
% Accelerations J2 [Equations 33]
% =========================================================================
if J2
    aR_J2 = 3 * mu * J2 * R^2 / (2 * B^8 * a^4 * G^2) * (12 * (Q1 * cos(L) - Q2 * sin(L))^2 - G^2) * Phi^4;
    aT_J2 = 12 * mu * J2 * R^2 / (B^8 * a^4 * G^2) * (Q2 * cos(L) + Q1 * sin(L)) * (Q1 * cos(L) - Q2 * sin(L)) * Phi^4;
    aN_J2 = 6 * mu * J2 * R^2 / (B^8 * a^4 * G^2) * (Q1 * cos(L) - Q2 * sin(L)) * (1 - Q1^2 - Q2^2) * Phi^4;
else
    aR_J2 = 0;
    aT_J2 = 0;
    aN_J2 = 0;
end

% =========================================================================
% Constant Acceleration in RTH [Equations 4]
% =========================================================================
if thrust.thrust_rth && (thrust.rth_flag_eps_r2 == 0)
    if ~eclipse
               
       % Acceleration - updated based on the mass
       thrust.epsilon_rth = thrust.thrust_rth / m;
        aR_rth = thrust.epsilon_rth * cos(thrust.beta_rth) * cos(thrust.alpha_rth);
        aT_rth = thrust.epsilon_rth * cos(thrust.beta_rth) * sin(thrust.alpha_rth);
        aN_rth = thrust.epsilon_rth * sin(thrust.beta_rth) ;
    elseif eclipse
        aR_rth = 0;
        aT_rth = 0;
        aN_rth = 0;
    end
else
    aR_rth = 0;
    aT_rth = 0;
    aN_rth = 0;
end


% =========================================================================
% Constant Acceleration in RTH [Equations 4]
% =========================================================================
if thrust.thrust_rth && (thrust.rth_flag_eps_r2 == 1)
    if ~eclipse
               
       % Acceleration - updated based on the mass
       thrust.epsilon_rth = thrust.thrust_rth / m;
        aR_rth_r2 = thrust.epsilon_rth/(r^2) * cos(thrust.beta_rth) * cos(thrust.alpha_rth);
        aT_rth_r2 = thrust.epsilon_rth/(r^2) * cos(thrust.beta_rth) * sin(thrust.alpha_rth);
        aN_rth_r2 = thrust.epsilon_rth/(r^2) * sin(thrust.beta_rth) ;
    elseif eclipse
        aR_rth_r2 = 0;
        aT_rth_r2 = 0;
        aN_rth_r2 = 0;
    end
else
    aR_rth_r2 = 0;
    aT_rth_r2 = 0;
    aN_rth_r2 = 0;
end



% =========================================================================
% Solar Radiation Pressure
% =========================================================================
if thrust.epsilon_inertial 
   if ~eclipse
       
       % Sun centered position of the Earth in the ecliptic reference frame
       [r_Sun_Earth, v_Earth] = EphSS(3, t0 + t * TU / 86400);
       r_Sun_Earth = r_Sun_Earth/DU;
       
       % Earth-Sun vector in the ecliptic reference frame
       r_Earth_Sun = - r_Sun_Earth;
       
       % Earth centered position of the Sun in the equatorial reference
       % frame: rotation from the ecliptic to the equatorial plane
       eps = astro_constants(8);
       r_Sun = [1 0 0; 0 cos(-eps) sin(-eps); 0 -sin(-eps) cos(-eps)] * r_Earth_Sun';
       
       % Spacecraft cartesian coordinates [DU]
       kep_SC = eq2kep([x(1:5); theta]);
       state_SC = kep2cart(kep_SC,mu);
       r_SC = state_SC(1:3);
       
       % SC/Sun distance [DU]
       r_SC_Sun = r_Sun' - r_SC;
       
       % Acceleration SRP in Earth centered frame []
        a_SRP = - P_Sun * Cr * A_m_SRP / m *  10^(-3) * ( TU^2 / DU) * r_SC_Sun / norm(r_SC_Sun)^3 * AU^2 / DU^2;
        
        % Acceleration SRP in rth reference frame
        aRTN_SRP = car_rthT(a_SRP, state_SC);
%         keyboard
    
        aR_SRP = aRTN_SRP(1);
        aT_SRP = aRTN_SRP(2);
        aN_SRP = aRTN_SRP(3);
        
%         aR_SRP = thrust.epsilon_inertial * cos(beta_Sun) * cos(alpha_Sun);
%         aT_SRP = thrust.epsilon_inertial * cos(beta_Sun) * sin(alpha_Sun);
%         aN_SRP = thrust.epsilon_inertial * sin(beta_Sun) ;
    elseif eclipse
        aR_SRP = 0;
        aT_SRP = 0;
        aN_SRP = 0;
    end
else
    aR_SRP = 0;
    aT_SRP = 0;
    aN_SRP = 0;
end

% =========================================================================
% Third body acceleration - Sun
% =========================================================================

if input.Sun
    
    % Transformation from MJD2000 to JD
    JD = mjd20002jd(t);

    % Sun position in Earth reference frame [km]
    r_Sun = sun_position(JD)';
    % Sun position in Earth reference frame [DU]
    r_Sun = r_Sun / DU;

    % Spacecraft cartesian coordinates [DU]
    kep_SC = eq2kep(x);
    state_SC = kep2cart(kep_SC,mu);
    r_SC = state_SC(1:3);

    % Sun acceleration in Earth reference frame [DU/TU^2]
    aSun = mu * ( (r_Sun - r_SC)/ norm(r_Sun - r_SC)^3 - r_Sun / norm(r_Sun)^3);

    % Transformation from Earth centered reference frame to RTN reference
    % frame
    aSun_rth = car_rthT(aSun,state_SC);
 
    aR_Sun = aSun_rth(1);
    aT_Sun = aSun_rth(2);
    aN_Sun = aSun_rth(3) ;
else
    aR_Sun = 0;
    aT_Sun = 0;
    aN_Sun = 0;
end

% =========================================================================
% Third body acceleration - Moon
% =========================================================================

if input.Moon
    
    % Transformation from MJD2000 to JD
    JD = mjd20002jd(t);

    % Sun position in Earth reference frame [km]
    r_Moon = moon_position(JD)';
    % Sun position in Earth reference frame [DU]
    r_Moon = r_Moon / DU;
    
    % Spacecraft cartesian coordinates [DU]
    kep_SC = eq2kep(x);
    state_SC = kep2cart(kep_SC,mu);
    r_SC = state_SC(1:3);

    % Sun acceleration in Earth reference frame [DU/TU^2]
    aMoon = mu * ( (r_Moon - r_SC)/ norm(r_Moon - r_SC)^3 - r_Moon / norm(r_Moon)^3);

    % Transformation from Earth centered reference frame to RTN reference
    % frame
    aMoon_rth = car_rthT(aMoon,state_SC);
 
    aR_Moon = aMoon_rth(1);
    aT_Moon = aMoon_rth(2);
    aN_Moon = aMoon_rth(3) ;
else
    aR_Moon = 0;
    aT_Moon = 0;
    aN_Moon = 0;
end

% =========================================================================
% Tangential Acceleration [Equations 25]
% =========================================================================
if thrust.thrust_t
    if ~eclipse
        % Acceleration - updated based on the mass
        thrust.epsilon_t = thrust.thrust_t / m;
        aR_Tang = thrust.epsilon_t * (P2 * sin(L) - P1 * cos(L))/D * cos(thrust.beta_t);
        aT_Tang = thrust.epsilon_t * (1 + P1 * sin(L) + P2 * cos(L))/D * cos(thrust.beta_t); 
        aN_Tang = thrust.epsilon_t * sin(thrust.beta_t);
    elseif eclipse
        aR_Tang = 0;
        aT_Tang = 0;
        aN_Tang = 0;
    end     
else
    aR_Tang = 0;
    aT_Tang = 0;
    aN_Tang = 0;
end

% =========================================================================
% Control Acceleration at Perigee and Apogee - follows the work of Zuiani
% =========================================================================
if thrust.flag_peri_apo
    
    dL_a    = control.dL_a;
    dL_p    = control.dL_p;
    beta_a  = control.beta_a;
    beta_p  = control.beta_p;
    ts      = control.ts;
    
    
    L_a_i   = interp1q_c(ts, dL_a, abs(t));
    L_p_i   = interp1q_c(ts, dL_p, abs(t));
    beta_a_i = interp1q_c(ts,beta_a,abs(t));
    beta_p_i = interp1q_c(ts,beta_p,abs(t));
    
    if L_a_i < 0
        k_a      = -k_a;
        L_a_i    = -L_a_i;
        beta_a_i = -beta_a_i;
    end
    
    if L_p_i < 0
        k_p     = -k_p;
        L_p_i   = -L_p_i;
        beta_p_i = -beta_p_i;
    end
    
    % Perigee thrust control
    epsilon_t = thrust.T_adim / m;
    aR_Tang_p = k_p * epsilon_t * (P2 * sin(L) - P1 * cos(L))/D * cos(beta_p_i);
    aT_Tang_p = k_p * epsilon_t * (1 + P1 * sin(L) + P2 * cos(L))/D * cos(beta_p_i);
    aN_Tang_p = k_p * epsilon_t * sin(beta_p_i);
    
    % Apogee thrust control
    epsilon_t = thrust.T_adim / m;
    aR_Tang_a = k_a * epsilon_t * (P2 * sin(L) - P1 * cos(L))/D * cos(beta_a_i);
    aT_Tang_a = k_a * epsilon_t * (1 + P1 * sin(L) + P2 * cos(L))/D * cos(beta_a_i);
    aN_Tang_a = k_a * epsilon_t * sin(beta_a_i);
    
    % Flag for perigee thrust arc. We are inside the pergee thrust arc if:
    % -> theta = L_p_i
    % -> theta = 360 - L_p_i
    % -> theta < L_p_i || theta > 360 - L_p_i
    perigee_flag = (sign(L_p_i -theta)==1 || ...
        sign(theta - 2*pi + L_p_i)==1 || ...
        (L_p_i  - theta == 0)  || (theta == 2*pi  - L_p_i));
    
    % Flag for the apogee thrust arc. We are in the apogee thrust arc is:
    % -> theta > 180 - L_a_i && theta < 180 + L_a_i
    apogee_flag  = (sign(theta-pi+L_a_i)==1 && ...
        sign(pi+L_a_i-theta)==1);
    
    aR_peri_apo = aR_Tang_p * perigee_flag + aR_Tang_a * apogee_flag;
    aT_peri_apo = aT_Tang_p * perigee_flag + aT_Tang_a* apogee_flag;
    aN_peri_apo = aN_Tang_p * perigee_flag + aN_Tang_a * apogee_flag;
    
    % Put this flag to 1 whenever the engine is on
    flag_thrust_peri_apo = perigee_flag || apogee_flag;
else
    flag_thrust_peri_apo = 0;
    aR_peri_apo = 0;
    aT_peri_apo = 0;
    aN_peri_apo = 0;
   
end
% =========================================================================
% Inertial Acceleration [Equations 21]
% =========================================================================
% aR_Inert = eps_Inert * cos(beta_Inert) * 

% =========================================================================
% Total accelerations
% =========================================================================
% keyboard
aR = aR_Drag + aR_J2 + aR_Tang + aR_rth + aR_rth_r2 + aR_SRP +...
    aR_Sun + aR_Moon + aR_peri_apo;
aT = aT_Drag + aT_J2 + aT_Tang + aT_rth +  aT_rth_r2 + aT_SRP +...
    aT_Sun + aT_Moon + aT_peri_apo;
aN = aN_Drag + aN_J2 + aN_Tang + aN_rth + aN_rth_r2 + aN_SRP + ...
    aN_Sun + aN_Moon + aN_peri_apo;

% keyboard
%% Equations of motion [Equations 2]

da_dt = (2/B) * sqrt(a^3/mu) * ( (P2 * sin(L) - P1 * cos(L)) * aR + Phi * aT);

dP1_dt = B * sqrt(a/mu) * ( - cos(L) * aR + ...
    ( ((P1 + sin(L))/Phi) + sin(L)) * aT - ...
    P2 * ((Q1 * cos(L) - Q2 * sin(L)) / Phi) * aN);
dP2_dt = B * sqrt(a/mu) * (   sin(L) * aR +...
    ( ((P2 + cos(L))/Phi) + cos(L)) * aT +...
    P1 * ((Q1 * cos(L) - Q2 * sin(L)) / Phi) * aN);

dQ1_dt = 0.5 * B * sqrt(a/mu) * (1 + Q1^2 + Q2^2) * sin(L) / Phi * aN;
dQ2_dt = 0.5 * B * sqrt(a/mu) * (1 + Q1^2 + Q2^2) * cos(L) / Phi * aN;
p_r = 1 + P1 * sin(L) + P2 * cos(L);
dl_dt = n - (r / h) * ( (a/(a+b) * (p_r) * (P1 * sin(L) + P2 * cos(L)) + (2*b/a)) * aR +...
                       a/(a+b) * (1+ (p_r)) * (P1 * cos(L) - P2 * sin(L)) * aT + ...
                       (Q1 * cos(L) - Q2 * sin(L)) * aN);
       
                   
dm_dt = -abs(thrust.thrust_t + thrust.thrust_rth + thrust.T_adim * flag_thrust_peri_apo) / (thrust.Isp * g0);                   
% dm_dt=0;
dEq_dt = [da_dt; dP1_dt; dP2_dt; dQ1_dt; dQ2_dt; dl_dt; dm_dt];

