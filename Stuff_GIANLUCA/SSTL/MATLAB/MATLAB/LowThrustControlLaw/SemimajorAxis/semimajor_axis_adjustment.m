function [Delta_V, T_Analytic, Kep_Analytic] = semimajor_axis_adjustment(initial_orbit, final_orbit, time_int, engine, m0, element_target, constants)

% Function for the variation of the semimajor axis using the optimal
% thrusting angles: azimuth directed as the flight path angle and elevation
% equal to zero. DeltaV can be computed by Edelbaum relationship for e~=0.
% For more accurate results the propagation can be stopped when the desired
% value of semimajor axis has been reached and the DeltaV computed
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

% In integrand_mod it will be assumed that there are no thurst arc at
% perigee and apogee but that during the entire rest of the orbit there is
% a constant tangential acceleration - therefore dL_p and dL_a are zero
k_a    = 0;
k_p    = 0;
alfa   = 0;
dL_a   = zeros(10,1);
dL_p   = zeros(10,1);
beta_a = zeros(10,1);
beta_p = zeros(10,1);
eta    = zeros(10,1);
csi    = zeros(10,1);
u_ratio = zeros(10,1);
ts =   linspace(0,ToF,10)';
t_wait = 0;
t_om_change = 0;
Dt_om_change = 0;
e_target = 0;
r_belt = 0;
T_Sun_adim = 0;
T_adim = 0;
t0 = 0;

% Thrust - tangential acceleration. Constant tangential acceleration rth is
% zero
thrust.thrust_t = T*1e-3 /constants.DU * constants.TU^2;
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
options_ODE = odeset('AbsTol',1e-12,'RelTol',1e-12,'events',@(t,x)events_stop_ode(t,x,element_target,final_orbit.a / constants.DU));
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
control.k_a    = k_a;
control.k_p    = k_p;
control.alfa   = alfa;
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
[T_Analytic,Equin_Analytic] = ode113(@(tt,xx)integrand_analytic_propagator(tt, xx, control, input, constants), ...
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
% DeltaV computed considering the final time [km/s]
Delta_V.DeltaV_time = T_Analytic(end)*constants.TU * f;

% Delta V from Edelbaum [km/s]
v_0 = sqrt(constants.mu_dim/initial_orbit.a);
v_f = sqrt(constants.mu_dim/final_orbit.a);
Delta_V.DeltaV_Edelbaum = sqrt(v_0^2 + v_f^2 - 2 * v_0 * v_f);

% DeltaV computed from rocket equation and mass consumption [m/s]
DeltaV = Isp * 9.8 * log(Equin_Analytic(1,7) / Equin_Analytic(end,7));
% DeltaV [km/s]
Delta_V.DeltaV = DeltaV * 1e-3;
