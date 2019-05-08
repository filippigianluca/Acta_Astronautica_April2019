% ========================================================================
% Main script for propagation of orbit in an Moon-centered reference system 
% using averaged analytical propagator. 
% Perturbations that can be included:
% - J2, J3, J4, J5
% - Solar radiation pressure (with eclipses)???
% - Third body (Sun)
% - Low-thrust tangential acceleration
% - Low-thrust constant acceleration in the rth reference frame
% =========================================================================
% Author: Marilena Di Carlo, 2016
% email: marilena.di-carlo@strath.ac.uk


clear
% close all
% clc

% Add folder to path (change this)
addpath(genpath('../spaceart_toolbox'))
addpath(genpath('IntegralsJ3'))



%% Constants

% Distance unit DU - Moon Radius [km]
constants.DU = 1737.5;

% mu [km^3/s^2]
mu_Moon = astro_constants(20);

% Time unit TU for Earth orbit [s]
constants.TU = sqrt(constants.DU^3 / mu_Moon);

%Moon Gravitational Constant [DU^3/TU^2]
constants.mu = 1;

% Sun gravitational constant [DU^3 / TU^2]
constants.mu_Sun = astro_constants(4) / constants.DU^3 * constants.TU^2;

% Standard free fall (the acceleration due to gravity on the
% Earth's surface at sea level) [m/s^2] 
constants.g0_dim = 9.8665;

% Standard free fall (the acceleration due to gravity on the
% Earth's surface at sea level) [DU/TU^2]
constants.g0 = constants.g0_dim * 1e-3 / constants.DU * constants.TU^2;

% Adimensional Earth radius [DU] (in this case is Moon)
constants.R_Earth = 1;

% J2, J3, J4, J5
constants.J2 =  202.43e-6;
constants.J3 =  8.476e-6;
constants.J4 =  0;
constants.J5 =  0;

% Seconds in a day
constants.sec_day = 86400;

% Fix the following:
% Solar radiation pressure [N/m^2 = kg /(m s2)] and [kg/(DU TU2)]
P_Sun = 4.56e-6;
constants.P_Sun = P_Sun;
% % In [kg /km s^2]
% P_Sun = P_Sun * 1000;
% constants.P_Sun = P_Sun * constants.DU * constants.TU^2;

% Radius of the radiation belt (from the original code of Federico)
r_belt = 0;


%% User input - MODIFY THIS

% Do you want to realise a comparison with numerical propagation? Possible
% only when some perturbations are used. 0 for no, 1 for yes
num_comparison = 0;

% Plot the results?
plot_flag = 1;

% -------------------------------------------------------------------------
% Initial date of the propagation
% -------------------------------------------------------------------------
date.year    = 2020;
date.month   = 3;
date.day     = 1;
date.hour    = 0;
date.minutes = 0;
date.seconds = 0;

% -------------------------------------------------------------------------
% Final date of the propagation
% -------------------------------------------------------------------------
date_end.year    = 2020;
date_end.month   = 3;
date_end.day     = 1;
date_end.hour    = 1;
date_end.minutes = 53;
date_end.seconds = 0;

% -------------------------------------------------------------------------
% Initial osculating keplerian elements (semimajor axis in DU, eccentricity,
% inclination, right ascension, perigee argument, true anomaly - angles in
% rad)
% -------------------------------------------------------------------------
% kep_osc0 = [(500+constants.DU)/ constants.DU ...
%              0 ...
%             30  * pi/180 ...
%             30 * pi/180 ...
%             30 * pi/180 ...
%               0 * pi/180];
          
%           
kep_osc0 = [(1788.6409)/ constants.DU ...
          6.3785e-4 ...
            59.9964  * pi/180 ...
            359.3139 * pi/180 ...
            180.1486 * pi/180 ...
              0 * pi/180];

% -------------------------------------------------------------------------
% Initial mass of the spacecraft [kg]
% -------------------------------------------------------------------------
m0 = 24;

% -------------------------------------------------------------------------
% Engine specific impulse [s]
% ------------------------------------------------------------------------
Isp = 750;

% -------------------------------------------------------------------------
% Engine Thrust [N]
% -------------------------------------------------------------------------
T =   5e-3;

% -------------------------------------------------------------------------
% Spacecraft characteristics drag and SRP parameters
% -------------------------------------------------------------------------
% Reflectivity coefficient
Cr = 1.3;

% Drag coefficient
CD = 2.2;

% Area to mass ratio for drag [km^2/kg]
A_m_DRAG = 1e-8;


% Area to mass ratio for Solar Radiation Pressure [km^2/kg]
A_m_SRP = 1e-8;


% -------------------------------------------------------------------------
% Perturbations
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Atmospheric drag?
flag_drag = 0;

% -------------------------------------------------------------------------
% J2, J3, J4, J5?
flag_J2 = 1;
flag_J3 = 1;
flag_J4 = 0;
flag_J5 = 0;

% -------------------------------------------------------------------------
% Earth flattening? (only for numerical propagation of the integrals of the
% drag)
flag_Earth_flat = 0;

% -------------------------------------------------------------------------
% 3rd body? Sun and Moon
flag_3rd_Sun  = 0;
flag_3rd_Moon = 0;

% -------------------------------------------------------------------------
% Constant acceleration in rth reference frame on the entire trajectory?
flag_rth  = 0;
% Elevation and azimuth angle of the acceleration [rad]
beta_rth  = pi/2;
alpha_rth = pi/4;

% -------------------------------------------------------------------------
% Constant tangential acceleration on the entire trajectory?
flag_t    = 0;
% Elevation angle of the tangential acceleration [rad] (remove this?)
beta_t    = 0*pi/4;

% -------------------------------------------------------------------------
% Acceleration only at perigee and apogee?
flag_peri_apo = 1;
% Number of nodes for the interpolation
n = 4;
% 1/0 for apogee and perigee raising/decrease
k_a    = flag_peri_apo*1;
k_p    = flag_peri_apo*1;
% Azimuth angle
alfa   = flag_peri_apo;
% Semiamplitude of the apogee or perigee thrust arc
dL_a   = flag_peri_apo * 00 * pi/180 * ones(n,1);
dL_p   = flag_peri_apo * 90 * pi/180 * ones(n,1);
% Elevation of the perigee and apogee thrust arc
beta_a = flag_peri_apo * 0 * ones(n,1);
beta_p = flag_peri_apo * 0 * ones(n,1);
% Shift of the perigee thrsut arc wrt perigee
% eta    = flag_peri_apo * 0*(-pi/2-kep_osc0(5)) * ones(n,1);
eta    = flag_peri_apo *240*pi/180 * ones(n,1);
% csi could also be 0, 0.5 or 1. The effect on the averaged propagation is
% on the choice of the initial point for the propagation over one orbit.
% Could be one single propagation starting from perigee, one single
% propagation starting from apogee or two split propagations one from
% perigee to apogee and the second from apogee to perigee
csi    = flag_peri_apo * 0.5 * ones(n,1);
%
u_ratio = zeros(n,1);

% -------------------------------------------------------------------------
% Constant inertial acceleration (SRP)?
flag_SRP  = 0;

% -------------------------------------------------------------------------
% Take into account eclipse for SRP? 
flag_ecl  = 0;


% -------------------------------------------------------------------------
% Options for solver inside integrator (not sure if this is still used)
% -------------------------------------------------------------------------
options_SOL = optimset('display','off','TolFun',1e-5);



% -------------------------------------------------------------------------
% Define when the integration should be stopped.
% When a certain eccentricity is reached or when a certain altitude is
% reached?
% -------------------------------------------------------------------------
stop_eccentricity = 0;
stop_altitude     = 1;

% -------------------------------------------------------------------------
% Define limit value for the eccentricity or the altitude (depending on the
% choice above) to stop the integration. Altitude in km
% -------------------------------------------------------------------------
% stop_eccentricity_value = 0.1;
stop_altitude_value     = 10;

% -------------------------------------------------------------------------
% Options for numerical integration
% -------------------------------------------------------------------------
options_ODE_NUM = odeset('AbsTol',1e-6,'RelTol',1e-6);

% Integration steps (number of step for the numerical integration and for
% the update of the integrand of the analytical averaging integration)
steps = 10000;


%% Initialisation

% -------------------------------------------------------------------------
% Initial time and propagation time
% ------------------------------------------------------------------------- 
% Initial modified julian date [MJD2000]
date.MJD2000 = date2mjd2000([date.year date.month date.day ...
    date.hour date.minutes date.seconds]);

% Initial time
t0 = date.MJD2000;

% Final modified julian date [MJD2000]
date_end.MJD2000 = date2mjd2000([date_end.year date_end.month date_end.day ...
    date_end.hour date_end.minutes date_end.seconds]);

% Time of flight
ToF = (date_end.MJD2000 - date.MJD2000 ) *  constants.sec_day / constants.TU;


% -------------------------------------------------------------------------
% Conversions
% -------------------------------------------------------------------------
% Area to mass ratio for drag [DU^2/kg]
A_m_DRAG = A_m_DRAG / constants.DU^2;

% Area to mass ratio [DU^2/kg]
A_m_SRP = A_m_SRP / constants.DU^2;

% Adimensional specific impulse
thrust.Isp = Isp / constants.TU;



% -------------------------------------------------------------------------
% Altitude at which spacecraft is considered re-entered and integration is
% stopped
% -------------------------------------------------------------------------
r_reentry = constants.R_Earth + stop_altitude_value/constants.DU;

% -------------------------------------------------------------------------
% Options for ODE averaged analytical integration
% -------------------------------------------------------------------------
if stop_altitude
    options_ODE = odeset('AbsTol',1e-7,'RelTol',1e-7,'events',@(t,x)events_reentry(t,x,r_reentry,'altitude'));
elseif stop_eccentricity
    options_ODE = odeset('AbsTol',1e-7,'RelTol',1e-7,'events',@(t,x)events_reentry(t,x,r_reentry,'eccentricity'));
end


% -------------------------------------------------------------------------
% Convert from osculating to mean elements (J2) and generate initial
% equinoctial mean elements
% -------------------------------------------------------------------------
if flag_J2
    
    % Convert osculating elements into mean elements.
    % The function osc3mean
    % requires the argument of perigee as 5th element
    kep_mean0_tmp = osc2mean([kep_osc0(1:3) kep_osc0(5) kep_osc0(4) kep_osc0(6)], constants.R_Earth, constants.J2);
    
    % Rearrange elements in kep_mean0 so that right ascension is in 4th
    % position and perigee argument is in 5h posiion
    kep_mean0    = kep_mean0_tmp;
    kep_mean0(4) = kep_mean0_tmp(5);
    kep_mean0(5) = kep_mean0_tmp(4);
   
%     kep_mean0(2)=0;
    % Initial mean equinoctial elements [DU, rad]
    Equin_start = kep2eq(kep_mean0);
 
   
else
    
    Equin_start = kep2eq(kep_osc0);
end

% -------------------------------------------------------------------------
% Create initial state vector for averaged analytical propagation
% -------------------------------------------------------------------------
% Initial vector for propagation (six equinoctial elements, mass, time below radiation
% belt, time in eclipse)
xx_0 = [Equin_start; m0; 0; 0];



% -------------------------------------------------------------------------
% Create initial state vector for numerical propagation
% -------------------------------------------------------------------------
if num_comparison
    % Initial osculating equinoctial elements [DU, rad]
    Equin_osc0    = kep2eq(kep_osc0);
    Equin_osc0(6) = mod(kep_osc0(4) + kep_osc0(5), 2*pi);
end




% Propellant mass flow rate [kg/s]
% m_rate = T / (Isp * constants.g0);

% !!!!!! Attention! In the original code of Federico m_rate is defined
% using the equation below:
m_rate = 1 / (Isp * constants.g0_dim);
% Look at integrand_analytic_propagator!!!!

% Adimensional propellant mass flow rate
m_rate = m_rate   * constants.DU /  constants.TU / 1e-3;



%% Chebyshev coefficients

% Coefficients for the interpolation of the atmospheric density (exponential model)
% using Chebyshev expansion of order 4.
% Depending on the considered altitude, different coefficients will be
% used. The coefficients have been obtained for different range of
% altitude, defined below

% Vector to define the height limits for each region in which the
% interpolation is realized
% h_crossing = [150 250 350 500 700 1000 2000 3000 4000 384000];
h_crossing = [110 125 150 250 350 500 700 1000 2000 3000 4000 384000];
h_crossing = h_crossing / constants.DU;

% Load coefficients
% load Chebyshev_coefficients.mat
load Chebyshev_coefficients_110km.mat




%% Perturbations and inputs to the analytical propagator

% =========================================================================
% Tangential thrust
% =========================================================================
if flag_t
    
    % Adimensional thrust [kg DU / TU^2]
    thrust.thrust_t = -T  * 10^(-3) * ( constants.TU^2  / constants.DU);
    
    % Initial mass [kg]
    thrust.m0 = m0;
    
    % Elevation angle
    thrust.beta_t = beta_t;
    
else
    thrust.epsilon_t = 0;
    thrust.beta_t    = 0;
    thrust.thrust_t  = 0;
end

% =========================================================================
% Constant Acceleration in rth reference frame
% =========================================================================
if flag_rth
    
    % Adimensional thrust [kg DU / TU^2]
    thrust.thrust_rth = T  * 10^(-3) * ( constants.TU^2  / constants.DU);
    
    % Initial mass [kg]
    thrust.m0 = m0;
    
    % Elevation angle [rad]
    thrust.beta_rth = beta_rth;
    
    % Azimuth angle [rad]
    thrust.alpha_rth = alpha_rth;
    
else
    thrust.epsilon_rth = 0;
    thrust.beta_rth = 0;
    thrust.alpha_rth = 0;
    thrust.thrust_rth = 0;
end


% =========================================================================
% Inertial Acceleration - check this 
% =========================================================================
if flag_SRP
    
    % SRP Acceleration [m/s^2] - modify to add the SC/SUn distance!
    a_inertial = constants.P_Sun * Cr * A_m_SRP ;
    
    % Adimensional acceleration [DU/TU^2]
    thrust.epsilon_inertial = a_inertial * 10^(-3) * ( constants.TU^2  / constants.DU);
    
        % [kg DU/TU^2]
    T_Sun_adim = m0 * thrust.epsilon_inertial;
%     thrust.epsilon_inertial = 0;
    

    
    % Elevation angle
    thrust.beta_In = 0;
    
    % Initial azimuth angle
    thrust.alpha_In = 0;
    
else
    
    thrust.epsilon_inertial = 0;
    thrust.beta_In = 0;
    thrust.alpha_In = 0;
    T_Sun_adim = 0;
    
end

% Collect all thrust input
input.thrust = thrust;

% Earth flattening
input.Earth_flattening = flag_Earth_flat;


% =========================================================================
% Thrust at perigee and apogee
% ========================================================================
if flag_peri_apo
    T_adim = T * 10^(-3) * ( constants.TU^2  / constants.DU);
    thrust.T_adim = T_adim;
    thrust.flag_peri_apo = 1;
else
    T_adim = 0;
    thrust.T_adim = T_adim;
    thrust.flag_peri_apo = 0;
end

% =========================================================================
% J2
% ========================================================================
if flag_J2
    input.geopotential.J2 = constants.J2;
else
    input.geopotential.J2 = 0;
end

% =========================================================================
% J3
% =========================================================================
if flag_J3
    input.geopotential.J3 = constants.J3;
else
    input.geopotential.J3 = 0;
end

% =========================================================================
% J4
% =========================================================================
if flag_J4
    input.geopotential.J4 = constants.J4;
else
    input.geopotential.J4 = 0;
end

% =========================================================================
% J5
% =========================================================================
if flag_J5
    input.geopotential.J5 = constants.J5;
else
    input.geopotential.J5 = 0;
end

% =========================================================================
% Drag
% =========================================================================
if flag_drag == 1
    
    % How to compute the integrals? 0 for numerical integration, 1 for
    % analytical integration. 0 if you are at low altitude, where there is
    % a coupling with J2
    drag.num_an     = 0;
    drag.CD         = CD;
    drag.A_m        = A_m_DRAG;
    drag.coeff      = coeff;
    drag.h_crossing = h_crossing;
    
else
    drag.CD = 0;
end

% Collect all the input for the atmospheric drag
input.drag   = drag;

% 
% Waiting time that Federico needed for some specific application
t_wait = 0;
% Other times that Federico needed for some specific application
t_om_change = 0;
Dt_om_change = 0;
% Targeted eccentricity - not applicable
e_target = 0;

% Input variables
input.T_adim       = T_adim;
input.m_rate       = m_rate;
input.T_Sun_adim   = T_Sun_adim;
input.r_belt       = r_belt;
input.Dt_om_change = Dt_om_change;
input.options_SOL  = options_SOL;
input.ecl_flag     = flag_ecl;
input.t_0          = t0;

input.flag_3rd_Sun  = flag_3rd_Sun;
input.flag_3rd_Moon = flag_3rd_Moon;

%% Low Thrust Control Parameters - 
% The propagation is realized using a function that can be used also for
% the optimization using a control problem in which two thrust arcs per
% orbits are considered - therefore several control parameters are defined.

% Time vector for interpolation
ts =   linspace(0,ToF,n)';

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


%% Averaged analytical integration

tic
[T_Analytic,Equin_Analytic] = ode113(@(t,x)integrand_analytic_propagator(t, x, control, input, constants), ...
                                    linspace(0, ToF, steps),...    % Propagation time
                                    xx_0,...            % Initial condition for the propagation
                                    options_ODE);       % Options propagation
toc

% DeltaV [m/s]
DeltaV = Isp * 9.8 * log(Equin_Analytic(1,7) / Equin_Analytic(end,7));






%% Numerical integration

if num_comparison
    
    % Settings for numerical propagation
    input_num_int.J2 = input.geopotential.J2;
    input_num_int.DU = constants.DU;
    input_num_int.TU = constants.TU;
    input_num_int.muadim = constants.mu;
    input_num_int.R = constants.R_Earth;
    input_num_int.Sun = flag_3rd_Sun;
    input_num_int.Moon = flag_3rd_Moon;
    input_num_int.t0 = t0;
    
    input_num_int.SRP.P_Sun = constants.P_Sun;
    input_num_int.SRP.Cr = Cr;
    input_num_int.SRP.A_m = A_m_SRP * 1e-6 / constants.DU^2;
    
    % Increase the step of the numerical integration by a factor steps_num
    steps_num = 1;
    
    % Numerical propagation
    % +====================================================================
    % Gauss equations with equinoctial elements
    % +====================================================================
    % Initial true longitue
    L0 = Equin_osc0(6);
    
    % Initial true anomaly
    theta0 = kep_osc0(6);
    
    % Initial eccentric anomaly
    sin_E_end = sqrt(1-kep_osc0(2)^2) * sin(theta0) / (1+kep_osc0(2)*cos(theta0));
    cos_E_end = (kep_osc0(2) + cos(theta0)) / (1+kep_osc0(2)*cos(theta0));
    E_end = atan2(sin_E_end, cos_E_end);
    
    % Initial mean anomaly
    M_end = E_end - kep_osc0(2) * sin(E_end);
    
    % Initial mean longitude
    l0 = mod(M_end + kep_osc0(4) + kep_osc0(5), 2*pi);
    
    
    %     options_ODE_NUM = odeset('AbsTol',1e-13,'RelTol',1e-13,'events',@(t,x)events_reentry(t,x,hp_reentry));
    tic
    
    % Integrate Gauss equations wrt time using 6th non approximate equation for l (mean longitude)
    [T_Numeric, Equin_Numeric] = ode113(@(t,x)equations_motion3(t, x, input_num_int, thrust, control, drag, constants.g0),...
        linspace(0, ToF, steps_num*steps),...
        [Equin_osc0(1:5); l0; m0],...
        options_ODE_NUM);
    
    toc
    
end

%% Results

% Initialization of keplerian and cartesian coordinates and element
Keplerian_Analytic = zeros(6, length(Equin_Analytic));
Cartesian_Analytic = zeros(6, length(Equin_Analytic));

% Compute keplerian and cartesina coordinates
for k = 1 : length(Equin_Analytic)
    
    % keplerian elements and cartesian coordinated from analytical solution
    Keplerian_Analytic(:,k) = eq2kep(Equin_Analytic(k,1:6));
    Cartesian_Analytic(:,k) = kep2cart(Keplerian_Analytic(:,k),constants.mu);
    
end

% Apogee and perigee hight from analytical solution
h_apogee_Analytic = Keplerian_Analytic(1,:) .* (1 + Keplerian_Analytic(2,:)) - constants.R_Earth;
h_perigee_Analytic = Keplerian_Analytic(1,:) .* (1 - Keplerian_Analytic(2,:)) - constants.R_Earth;


if plot_flag
    
    % -------------------------------------------------------------------------
    % Apogee and perigee
    % ------------------------------------------------------------------------- 
    ff=figure(1);
    plot((T_Analytic(:)-T_Analytic(1)) * constants.TU /constants.sec_day, h_perigee_Analytic * constants.DU,'b','LineWidth',2)
    hold on
    plot((T_Analytic(:)-T_Analytic(1)) * constants.TU/constants.sec_day, h_apogee_Analytic*constants.DU,'r','LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('Height [km]')
    legend('Perigee height','Apogee height')
%     title('Analytic')
    
    
    % -------------------------------------------------------------------------
    % Keplerian elements 
    % -------------------------------------------------------------------------
    
    % Eccentricity and inclination
    e_Analytic = sqrt(Equin_Analytic(:,2).^2 + Equin_Analytic(:,3).^2);
    incl_Analytic = 2 * atan(sqrt(Equin_Analytic(:,4).^2 + Equin_Analytic(:,5).^2));
    
    
    
    figure('units','normalized','outerposition',[0 0 0.6 0.6])
    subplot(2,2,1)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/constants.sec_day, Keplerian_Analytic(1,:) * constants.DU,'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('Semimajor axis [km]')
    subplot(2,2,2)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/constants.sec_day, Keplerian_Analytic(2,:),'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('Eccentricity')
    subplot(2,3,4)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/constants.sec_day, Keplerian_Analytic(3,:)*180/pi,'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('Inclination [deg]')
    subplot(2,3,6)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/constants.sec_day, Keplerian_Analytic(5,:)*180/pi,'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('\omega [deg]')
    subplot(2,3,5)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/constants.sec_day, Keplerian_Analytic(4,:)*180/pi,'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('\Omega [deg]')
    set(gcf,'NextPlot','add');
    axes;
%     h = title('Analytical Propagation');
    set(gca,'Visible','off');
%     set(h,'Visible','on');
    
    
    % -------------------------------------------------------------------------
    % Equinoctial elements 
    % -------------------------------------------------------------------------
    figure('units','normalized','outerposition',[0 0 0.6 0.6])
    subplot(2,2,1)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/constants.sec_day, Equin_Analytic(:,1) * constants.DU,'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('Semimajor axis [km]')
    subplot(2,2,2)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/constants.sec_day, Equin_Analytic(:,2),'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('P_1')
    subplot(2,3,4)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/constants.sec_day, Equin_Analytic(:,3),'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('P_2')
    subplot(2,3,6)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/constants.sec_day, Equin_Analytic(:,4)*180/pi,'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('Q_1 [deg]')
    subplot(2,3,5)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/constants.sec_day, Equin_Analytic(:,5)*180/pi,'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('Q_2')
    
    
    
    
    % -------------------------------------------------------------------------
    % 3D plot
    % -------------------------------------------------------------------------
%     figure
%     % Initial trajectory
%     plot_trajectory(Cartesian_Analytic(1:3,1)*constants.DU,...
%         Cartesian_Analytic(4:6,1)*constants.DU/constants.TU,...
%         2*pi*sqrt((Keplerian_Analytic(1,1))^3/constants.mu)*constants.TU,...
%         constants.mu*constants.DU^3/constants.TU^2,'k',4);
%     hold on
%     % Final trajectory
%     plot_trajectory(Cartesian_Analytic(1:3,end)*constants.DU,...
%         Cartesian_Analytic(4:6,end)*constants.DU/constants.TU,...
%         2*pi*sqrt((Keplerian_Analytic(1,end))^3/constants.mu)*constants.TU,...
%         constants.mu*constants.DU^3/constants.TU^2,'b',4);
%     hold on
%     grid on
%     quiver3(0,0,0,5*constants.DU,0,0,'Color','k')
%     quiver3(0,0,0,0,5*constants.DU,0,'Color','k')
%     quiver3(0,0,0,0,0,5*constants.DU,'Color','k')
%       
%     for i = 2 : 50 : length(Keplerian_Analytic) -1
%         
%         plot_trajectory(Cartesian_Analytic(1:3,i)*constants.DU,...
%             Cartesian_Analytic(4:6,i)*constants.DU/constants.TU, ...
%             2*pi*sqrt((Keplerian_Analytic(1,i))^3/constants.mu)*constants.TU,...
%             constants.mu*constants.DU^3/constants.TU^2,[0.5 0.5 0.5],1);
%     end
%     hold on
%     [x_s,y_s,z_s] = sphere();
%     C = ones(size(z_s))*0;
%     hSurface = surf( x_s*constants.DU, y_s*constants.DU,z_s*constants.DU ,C) ;
%     set(hSurface,'FaceColor',[0 0 0],'FaceAlpha',0);
%     grid on
%     xlabel('x [km]','FontSize',12)
%     ylabel('y [km]','FontSize',12)
%     zlabel('z [km]','FontSize',12)
%     legend('Initial Orbit','Final Orbit')
%     axhand=gca;
%     set(axhand,'FontName','Helvetica','FontSize',12)
%     view(3)
%     axis equal

    
end


if num_comparison
    
    % Difference between analytic and numerical solution vector
    delta = Equin_Analytic(end,1:6) - Equin_Numeric(end,1:6);
    
    Keplerian_Numeric = zeros(6, steps);
    Cartesian_Numeric = zeros(6, steps);
    
    for k = 1 : steps
        
        
        % Keplerian elements and cartesian coordinates from numerical solution
        Keplerian_Numeric(:,k) = eq2kep(Equin_Numeric(k*steps_num,1:6));
        Cartesian_Numeric(:,k) = kep2cart(Keplerian_Numeric(:,k),constants.mu);
        
    end
    
    % Apogee and perigee hight from numerical solution
    h_apogee_Numeric = Keplerian_Numeric(1,:) .* (1 + Keplerian_Numeric(2,:)) - constants.R_Earth;
    h_perigee_Numeric = Keplerian_Numeric(1,:) .* (1 - Keplerian_Numeric(2,:)) - constants.R_Earth;
    
    
    figure
    subplot(1,2,1)
    plot((T_Numeric(:)-T_Numeric(1)) * constants.TU / 86400, h_perigee_Numeric * constants.DU,'r')
    hold on
    plot((T_Numeric(:)-T_Numeric(1))*constants.TU/constants.sec_day,   h_apogee_Numeric*constants.DU,'r')
    grid on
    xlabel('Time [days]')
    ylabel('Height [km]')
    legend('Perigee height','Apogee height')
    title('Numeric')
    

    
    
    %% Plot of the results
    
    % -------------------------------------------------------------------------
    % Apogee and perigee variations of analytically averaged and numerical
    % integration - plotted versus time elapsed from initial time
    % -------------------------------------------------------------------------
    figure
    subplot(1,2,1)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU / constants.sec_day, h_perigee_Numeric * constants.DU,'r')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU / constants.sec_day, h_perigee_Analytic * constants.DU,'b','LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('Perigee height [km]')
    legend('Numerical','Averaged analytical')
    
    subplot(1,2,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/constants.sec_day,   h_apogee_Numeric*constants.DU,'r')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/constants.sec_day, h_apogee_Analytic*constants.DU,'b','LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('Apogee height [km]')
    legend('Numerical','Averaged analytical')
    
    
    
    % -------------------------------------------------------------------------
    % Equinoctial elements for numerical and averaged analytical solution: a,
    % P1, P2 and Q1
    % -------------------------------------------------------------------------
    figure
    title('Equinoctial elements')
    subplot(2,2,1)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/constants.sec_day, Equin_Analytic(:,1)*constants.DU,'b')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU / constants.sec_day,Equin_Numeric(1:steps_num:end,1)*constants.DU,'r')
    xlabel('Time [days]')
    ylabel('a [km]')
    legend('Averaged analytical','Numerical')
    
    subplot(2,2,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU / constants.sec_day, Equin_Analytic(:,2),'b')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/constants.sec_day, Equin_Numeric(1:steps_num:end,2),'r')
    xlabel('Time [days]')
    ylabel('P_1')
    legend('Averaged analytical','Numerical')
    
    subplot(2,2,3)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Equin_Analytic(:,3),'b')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Equin_Numeric(1:steps_num:end,3),'r')
    xlabel('Time [days]')
    ylabel('P_2')
    legend('Averaged analytical','Numerical')
    
    subplot(2,2,4)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Equin_Analytic(:,4),'b')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Equin_Numeric(1:steps_num:end,4),'r')
    xlabel('Time [days]')
    ylabel('Q_1')
    legend('Averaged analytical','Numerical')
    
    
    figure
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Equin_Analytic(:,5),'b')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Equin_Numeric(1:steps_num:end,5),'r')
    xlabel('Time [days]')
    ylabel('Q_2')
    legend('Averaged analytical','Numerical')
    
    
    % -------------------------------------------------------------------------
    % Keplerian elements for numerical and averaged analytical solution: a,
    % e, i, RAAN, omega
    % -------------------------------------------------------------------------
    
    % Eccentricity analytical and numericla solution
    e_Numeric  = sqrt(Equin_Numeric(1:steps_num:end,2).^2  + Equin_Numeric(1:steps_num:end,3).^2);
    
    incl_Numeric  = 2 * atan (sqrt(Equin_Numeric(1:steps_num:end,4).^2 + Equin_Numeric(1:steps_num:end,5).^2));
    
    
    % -------------------------------------------------------------------------
    % Semimajor axis and eccentircity
    figure
    title('Keplerian elements')
    subplot(2,1,1)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU / constants.sec_day,Equin_Numeric(1:steps_num:end,1)*constants.DU,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/constants.sec_day, Equin_Analytic(:,1)*constants.DU,'b','LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('a [km]')
    legend('Numerical','Averaged analytical')
    
    subplot(2,1,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/constants.sec_day, e_Numeric,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU /constants.sec_day, e_Analytic,'b','LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('e')
    legend('Numerical','Averaged analytical')
    
    % -------------------------------------------------------------------------
    % Inclination, RAAN and perigee argument
    figure
    subplot(2,2,1)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, incl_Numeric*180/pi,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, incl_Analytic*180/pi,'b','LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('i [deg]')
    legend('Numerical','Averaged analytical')
    
    subplot(2,2,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Keplerian_Numeric(4,:)*180/pi,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Keplerian_Analytic(4,:)*180/pi,'b','LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('\Omega [deg]')
    legend('Numerical','Averaged analytical')
    
    subplot(2,2,3)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Keplerian_Numeric(5,:)*180/pi,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Keplerian_Analytic(5,:)*180/pi,'b','LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('\omega [deg]')
    legend('Numerical','Averaged analytical')
    
    
    
        figure('units','normalized','outerposition',[0 0 0.6 0.6])
    subplot(2,2,1)
       plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU / constants.sec_day,Equin_Numeric(1:steps_num:end,1)*constants.DU,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/constants.sec_day, Equin_Analytic(:,1)*constants.DU,'b','LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('a [km]')
    legend('Numerical Propagation','Averaged Analytical Propagation')
    axhand=gca;
set(axhand,'FontSize',14)
    subplot(2,2,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/constants.sec_day, e_Numeric,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU /constants.sec_day, e_Analytic,'b','LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('e')
        axhand=gca;
set(axhand,'FontSize',14)
%     legend('Numerical','Averaged analytical')
    subplot(2,3,4)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, incl_Numeric*180/pi,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, incl_Analytic*180/pi,'b','LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('i [deg]')
        axhand=gca;
set(axhand,'FontSize',14)
%     legend('Numerical','Averaged analytical')
    subplot(2,3,6)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Keplerian_Numeric(4,:)*180/pi,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Keplerian_Analytic(4,:)*180/pi,'b','LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('\Omega [deg]')
        axhand=gca;
set(axhand,'FontSize',14)
%     legend('Numerical','Averaged analytical')
    subplot(2,3,5)
 plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, (Keplerian_Analytic(5,:)-Keplerian_Numeric(5,:))*180/pi,'k','LineWidth',2)
plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Keplerian_Numeric(5,:)*180/pi,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Keplerian_Analytic(5,:)*180/pi,'b','LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('\omega [deg]')
        axhand=gca;
set(axhand,'FontSize',14)
%     legend('Numerical','Averaged analytical')
    
    
    
    
    % -------------------------------------------------------------------------
    % DELTA of the Keplerian elements for numerical and averaged analytical solution:
    % a, e, i, RAAN, omega
    % -------------------------------------------------------------------------
    figure
    subplot(2,1,1)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/constants.sec_day, (Equin_Analytic(:,1)-Equin_Numeric(1:steps_num:end,1))*constants.DU,'k','LineWidth',2)
    xlabel('Time [days]')
    ylabel('\Delta a [km]')
    subplot(2,1,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU /constants.sec_day, e_Analytic-e_Numeric,'k','LineWidth',2)
    xlabel('Time [days]')
    ylabel('\Delta e')
    
    
    figure
    subplot(2,2,1)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, (incl_Analytic-incl_Numeric)*180/pi,'k','LineWidth',2)
    xlabel('Time [days]')
    ylabel('\Delta i [deg]')
    subplot(2,2,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, (Keplerian_Analytic(4,:)-Keplerian_Numeric(4,:))*180/pi,'k','LineWidth',2)
    xlabel('Time [days]')
    ylabel('\Delta \Omega [deg]')
    subplot(2,2,3)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, (Keplerian_Analytic(5,:)-Keplerian_Numeric(5,:))*180/pi,'k','LineWidth',2)
    xlabel('Time [days]')
    ylabel('\Delta \omega [deg]')
    
    
    
    
    
    figure('units','normalized','outerposition',[0 0 0.6 0.6])
    subplot(2,2,1)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/constants.sec_day, (Equin_Analytic(:,1)-Equin_Numeric(1:steps_num:end,1))*constants.DU,'k','LineWidth',2)
    xlabel('Time [days]')
    ylabel('\Delta a [km]')
    subplot(2,2,2)
   plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU /constants.sec_day, e_Analytic-e_Numeric,'k','LineWidth',2)
    xlabel('Time [days]')
    ylabel('\Delta e')
    subplot(2,3,4)
   plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, (incl_Analytic-incl_Numeric)*180/pi,'k','LineWidth',2)
    xlabel('Time [days]')
    ylabel('\Delta i [deg]')
    subplot(2,3,6)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, (Keplerian_Analytic(4,:)-Keplerian_Numeric(4,:))*180/pi,'k','LineWidth',2)
    xlabel('Time [days]')
    ylabel('\Delta \Omega [deg]')
    subplot(2,3,5)
 plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, (Keplerian_Analytic(5,:)-Keplerian_Numeric(5,:))*180/pi,'k','LineWidth',2)
    xlabel('Time [days]')
    ylabel('\Delta \omega [deg]')
    
    
    
    
    
    
end
% Keplerian_Analytic(1,end)*(1-Keplerian_Analytic(2,end))*constants.DU-constants.DU
% Keplerian_Analytic(1,end)*(1+Keplerian_Analytic(2,end))*constants.DU-constants.DU

% dV(iijk) = DeltaV;
%  time(iijk) = (T_Analytic(end)-T_Analytic(1))*constants.TU/86400;
%  mass(iijk) = - (Equin_Analytic(end,7)-Equin_Analytic(1,7));
 


