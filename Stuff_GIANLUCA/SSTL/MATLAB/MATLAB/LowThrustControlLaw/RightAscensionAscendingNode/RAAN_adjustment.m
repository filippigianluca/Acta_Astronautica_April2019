function [Delta_V, T_Analytic, Kep_Analytic] = RAAN_adjustment(initial_orbit, final_orbit, time_int, engine, m0, element_target, k_p, k_a, constants)

% Function for the variation of the RAAN using the optimal
% thrusting angles: azimuth zero and elevation pi/2 changing sign twice per orbit. 
% The propagation stops when the desired
% value of inclination has been reached and the DeltaV is computed
% accordingly
% Input: initial_orbit -> structure with the elements of the initial orbit
%        final_orbit   -> structure with the elements of the final orbit
%        time_int      -> maximum integration time
%        engine        -> structure with the properties of the electric engine
%        m0            -> initial mass of the spacecraft
%        element_target-> string with the definition of the element that is
%                         adjusted
%        constants     ->

% Marilena Di Carlo, 2015
% marilena.di-carlo@strath.ac.uk

% Engine: thrust [kg DU/TU^2], specific impulse [TU] and acceleration [DU/TU^2]
T   = engine.T;
Isp = engine.Isp;
f   = engine.f;

% Time of integration
ToF    = time_int / constants.TU;

% In integrand_mod it will be assumed that there are thurst arcs at perigee
% and apogee of semiamplitude pi/2 and that the center of this thrust arc
% can be shifted by omega to account for sign(sin(omega+theta))
% k_a    = -1;
% k_p    = -1;
alfa   = 0;
n = 10;
% Semiamplitude thrust arcs
dL_a   = pi/2*ones(n,1);
dL_p   = pi/2*ones(n,1);
beta_a = -k_a * pi/2*ones(n,1);
beta_p = k_p * pi/2*ones(n,1);
eta    = (pi/2-initial_orbit.omega) * ones(n,1);
% keyboard
csi    = zeros(n,1);
u_ratio = zeros(n,1);
ts =   linspace(0,ToF,n)';
t_wait = 0;
t_om_change = 0;
Dt_om_change = 0;
e_target = 0;
r_belt = 0;
T_Sun_adim = 0;
T_adim = T * 1e-3 / constants.DU  * constants.TU^2;
t0 = 0;

% Thrust - tangential acceleration. Constant tangential acceleration rth is
% zero
thrust.thrust_t = 0;
thrust.beta_rth = 0;
thrust.alpha_rth = 0;
thrust.thrust_rth = 0;
thrust.Isp = Isp/constants.TU ;

% Initial mass [kg]
thrust.m0 = m0;

% Elevation angle - zero
thrust.beta_t = 0;

% Collect all thrust input
input.thrust = thrust;

% No drag
drag.CD = 0;

% Adimensional mass rate
% m_rate = thrust.thrust_t / ( thrust.Isp * constants.g0);
m_rate = 1 / (Isp * constants.g0_m_s);
m_rate = m_rate   * constants.DU /  constants.TU / 1e-3;

flag_ecl = 0;
steps = 3000;
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,'events',@(t,x)events_stop_ode_eq(t,x,element_target,final_orbit.Omega));
options_SOL = optimset('display','off','TolFun',1e-5);

% Initial conditions in keplerian coordinates
xx_0 = [initial_orbit.a/constants.DU initial_orbit.e initial_orbit.i initial_orbit.Omega initial_orbit.omega 0];

% Initial conditions in equinoctial coordinates
xx_0 = kep2eq(xx_0);
xx_0 = [xx_0; m0; 0; 0];

% -------------------------------------------------------------------------
% Averaged analytical integration
% -------------------------------------------------------------------------

% ------- Control
% Time vector for interpolation
% ts =   linspace(0,ToF,n)';

% Control parameters
control.ts     = ts;
control.dL_a   = dL_a;
control.dL_p   = dL_p;
control.k_a    = 0;
control.k_p    = 0;
control.alfa_a = 0*ones(n,1);  
control.alfa_p = 0*ones(n,1);
control.beta_a = beta_a;
control.beta_p = beta_p;
control.eta    = eta;
control.csi    = csi;
control.u_ratio= u_ratio;
control.ts     = ts;


% ------ Input
input.geopotential.J2 = 0;
input.geopotential.J3 = 0;
input.geopotential.J4 = 0;
input.geopotential.J5 = 0;
input.drag.CD         = 0;
input.T_adim       = T_adim;
input.m_rate       = m_rate;
input.T_Sun_adim   = T_Sun_adim;
input.r_belt       = r_belt;
input.Dt_om_change = Dt_om_change;
input.options_SOL  = options_SOL;
input.ecl_flag     = flag_ecl;
input.t_0          = t0;
input.flag_3rd_Sun  = 0;
input.flag_3rd_Moon = 0;
input.Earth_flattening = 0;

% % Averaged analytical integration
% [T_Analytic,Equin_Analytic] = ode113(@(tt,xx) integrand_mod(tt, xx, ts, ...
%                         dL_a, dL_p, k_a, k_p, alfa, beta_a, beta_p, eta, csi, u_ratio, ...
%                         1, constants.J2, 1, constants.g0_DUTU, drag, thrust,...
%                         T_adim, m_rate, ...
%                         T_Sun_adim, e_target, r_belt, ...
%                         t0 + t_wait, t0 + t_wait + t_om_change, Dt_om_change, options_SOL, flag_ecl),...
%                         linspace(0, ToF, steps),...    % Propagation time
%                         xx_0,...            % Initial condition for the propagation
%                         options_ODE);       % Options propagation
%      keyboard             
[T_Analytic,Equin_Analytic] = ode113(@(tt,xx)integrand_analytic_propagator_new(tt, xx, control, input, constants), ...
    linspace(0, ToF, steps),...    % Propagation time
    xx_0,...            % Initial condition for the propagation
    options_ODE);       % Options propagation);


                    
                    
% Transformation from equinoctial to keplerian elements
for i = 1 : length(Equin_Analytic)
    Kep_Analytic(i,:) = eq2kep(Equin_Analytic(i,1:6));
end


% =========================================================================
% DeltaV
% =========================================================================
Delta_V = Isp * 9.8 * log(Equin_Analytic(1,7) / Equin_Analytic(end,7))*1e-3;
