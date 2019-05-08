% ========================================================================
% Main script for propagation of orbit in an Sun-centered reference system 
% using averaged analytical propagator. 
% Perturbations that can be included:
% - Low-thrust tangential acceleration
% - Low-thrust constant acceleration in the rth reference frame
% =========================================================================
% Author: Marilena Di Carlo, 2016
% email: marilena.di-carlo@strath.ac.uk


clear
close all
clc

% Add folder to path (change this)
addpath(genpath('../spaceart_toolbox'))



%% Constants

% Time unit TU for interplanetary orbit [days]
constants.TU = 58.13;

% Distance unit DU - Earth-Sun distance [km]
constants.DU = astro_constants(2);

% Sun Gravitational Constant [DU^3/TU^2]
constants.mu = 1;

% Seconds in a day
constants.sec_day = 86400;

% Standard free fall (the acceleration due to gravity on the
% Earth's surface at sea level) [m/s^2] 
constants.g0_dim = 9.8665;

% Standard free fall (the acceleration due to gravity on the
% Earth's surface at sea level) [DU/TU^2]
constants.g0 = constants.g0_dim * 1e-3 / constants.DU * (constants.TU * constants.sec_day)^2;

% Adimensional Earth radius [DU]
constants.R_Earth = 6378.136 / constants.DU;



% Radius of the radiation belt (from the original code of Federico)
r_belt = 0;

% Fix the following:
% Solar radiation pressure [N/m^2 = kg /(m s2)] and [kg/(DU TU2)]
P_Sun = 4.56e-6;
constants.P_Sun = P_Sun;
% % In [kg /km s^2]
% P_Sun = P_Sun * 1000;
% constants.P_Sun = P_Sun * constants.DU * constants.TU^2;



%% User input - MODIFY THIS

% Do you want to realise a comparison with numerical propagation? Possible
% only when some perturbations are used. 0 for no, 1 for yes
num_comparison = 1;

% Plot the results?
plot_flag = 1;

% -------------------------------------------------------------------------
% Initial date of the propagation
% -------------------------------------------------------------------------
date.year    = 2030;
date.month   = 3;
date.day     = 21;
date.hour    = 0;
date.minutes = 0;
date.seconds = 0;

% -------------------------------------------------------------------------
% Final date of the propagation
% -------------------------------------------------------------------------
date_end.year    = 2030;
date_end.month   = 12;
date_end.day     = 21;
date_end.hour    = 0;
date_end.minutes = 0;
date_end.seconds = 0;

% -------------------------------------------------------------------------
% Initial osculating keplerian elements (semimajor axis in DU, eccentricity,
% inclination, right ascension, perigee argument, true anomaly - angles in
% rad)
% -------------------------------------------------------------------------
kep_osc0 = [ 1.5116 ...
             0 ...
              80 * pi/180 ...
              10 * pi/180 ...
            60 * pi/180 ...
             200 * pi/180];

% -------------------------------------------------------------------------
% Initial mass of the spacecraft [kg]
% -------------------------------------------------------------------------
m0 = 2000;

% -------------------------------------------------------------------------
% Engine specific impulse [s]
% ------------------------------------------------------------------------
Isp = 4000;

% -------------------------------------------------------------------------
% Engine Thrust [N]
% -------------------------------------------------------------------------
T =   0.1;

% -------------------------------------------------------------------------
% Spacecraft characteristics drag and SRP parameters
% -------------------------------------------------------------------------
% Reflectivity coefficient
Cr = 1.3;

% Area to mass ratio for Solar Radiation Pressure [km^2/kg]
A_m_SRP = 1e-8;


% -------------------------------------------------------------------------
% Perturbations
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Constant acceleration in rth reference frame on the entire trajectory?
flag_rth  = 0;
% Elevation and azimuth angle of the acceleration [rad]
beta_rth  = 10*pi/180;
alpha_rth = pi/2;

% -------------------------------------------------------------------------
% Constant acceleration in rth reference frame on the entire trajectory,
% with magnitude as 1/r^2?
flag_rth_r2  = 1;
% Elevation and azimuth angle of the acceleration [rad]
beta_rth_r2  = 0*15*pi/180;
alpha_rth_r2 = pi/2;
R_ref = 1;

% -------------------------------------------------------------------------
% Constant tangential acceleration on the entire trajectory?
flag_t    = 0 ;
% Elevation angle of the tangential acceleration [rad]
beta_t    = 0;

% -------------------------------------------------------------------------
% Acceleration only at perigee and apogee?
flag_peri_apo = 0;
% Number of nodes for the interpolation
n = 4;
% 1/0 for apogee and perigee raising/decrease
k_a    = 0;
k_p    = 0;
% Azimuth angle
alfa   = 0;
% Semiamplitude of the apogee or perigee thrust arc
dL_a   = flag_peri_apo*90*pi/180 * ones(n,1);
dL_p   = flag_peri_apo*90*pi/180 * ones(n,1);
% Elevation of the perigee and apogee thrust arc
beta_a = flag_peri_apo*pi/2 * ones(n,1);
beta_p = flag_peri_apo*pi/2 * ones(n,1);
% Shift of the perigee thrsut arc wrt perigee
eta    = 0 * ones(n,1);
% csi could also be 0, 0.5 or 1. The effect on the averaged propagation is
% on the choice of the initial point for the propagation over one orbit.
% Could be one single propagation starting from perigee, one single
% propagation starting from apogee or two split propagations one from
% perigee to apogee and the second from apogee to perigee
csi    = flag_peri_apo * 0.5 * ones(n,1);
% csi = 0.5* ones(n,1);
%
u_ratio = zeros(n,1);

% -------------------------------------------------------------------------
% Constant inertial acceleration (SRP)?
flag_SRP  = 0;

% -------------------------------------------------------------------------
% Take into account eclipse for SRP? 
flag_ecl  = 0;

% -------------------------------------------------------------------------
% 3rd body? Moon
flag_3rd_Moon = 0;
flag_3rd_Sun  = 0;

% -------------------------------------------------------------------------
% Options for solver inside integrator (not sure if this is still used)
% -------------------------------------------------------------------------
options_SOL = optimset('display','off','TolFun',1e-5);



% -------------------------------------------------------------------------
%Stop integration when main belt is reached
% -------------------------------------------------------------------------
semimajor_main_belt = 2.85;
eccentricity_main_belt = 0.009;
apogee = 8;

% -------------------------------------------------------------------------
% Options for ODE averaged analytical integration
% -------------------------------------------------------------------------
% options_ODE = odeset('AbsTol',1e-7,'RelTol',1e-7,'events',@(t,x)events_main_belt(t,x,semimajor_main_belt,'semimajoraxis'));
% options_ODE = odeset('AbsTol',1e-7,'RelTol',1e-7,'events',@(t,x)events_main_belt(t,x,eccentricity_main_belt,'eccentricity'));
options_ODE = odeset('AbsTol',1e-7,'RelTol',1e-7,'events',@(t,x)events_main_belt(t,x,apogee,'apogee'));

% options_ODE   = odeset('AbsTol',1e-7,'RelTol',1e-7);
% -------------------------------------------------------------------------
% Options for numerical integration
% -------------------------------------------------------------------------
options_ODE_NUM = odeset('AbsTol',1e-13,'RelTol',1e-13);

% Integration steps (number of step for the numerical integration and for
% the update of the integrand of the analytical averaging integration)
steps = 50000;


%% Initialisation

% -------------------------------------------------------------------------
% Earth flattening? (only for numerical propagation of the integrals of the
% drag) - 0 for interplanetary transfers
flag_Earth_flat = 0;
% Earth flattening
input.Earth_flattening = flag_Earth_flat;

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
ToF = (date_end.MJD2000 - date.MJD2000 ) / constants.TU;


% -------------------------------------------------------------------------
% Conversions
% -------------------------------------------------------------------------
% Area to mass ratio [DU^2/kg]
A_m_SRP = A_m_SRP / constants.DU^2;

% Adimensional specific impulse
thrust.Isp = Isp / (constants.TU * constants.sec_day);



% -------------------------------------------------------------------------
% Generate initial equinoctial mean elements
% -------------------------------------------------------------------------
Equin_start = kep2eq(kep_osc0);


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
m_rate = m_rate   * constants.DU / (constants.TU * constants.sec_day) / 1e-3;




%% Perturbations and inputs to the analytical propagator

% =========================================================================
% Tangential thrust
% =========================================================================
if flag_t
    
    % Adimensional thrust [kg DU / TU^2]
    thrust.thrust_t = T  * 10^(-3) * (  (constants.TU*constants.sec_day)^2  / constants.DU);
    
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
    thrust.thrust_rth = T  * 10^(-3) * (  (constants.TU*constants.sec_day)^2  / constants.DU);
    
    thrust.rth_flag_eps_r2 = 0;
    
    % Initial mass [kg]
    thrust.m0 = m0;
    
    % Elevation angle [rad]
    thrust.beta_rth = beta_rth;
    
    % Azimuth angle [rad]
    thrust.alpha_rth = alpha_rth;
    

elseif flag_rth_r2
    
    % Adimensional thrust [kg DU / TU^2]
    % Here the thrust must be T * R_ref^2, where R_ref is the distance from
    % the Sun where T is equal to its maximum value. If DU=AU and the
    % maximum distance is at 1 AU, then the thrust is simply T. 
    thrust.thrust_rth = T  * 10^(-3) * (  (constants.TU*constants.sec_day)^2  / constants.DU);
    
    thrust.rth_flag_eps_r2 = 1;
    
    % Initial mass [kg]
    thrust.m0 = m0;
    
    % Elevation angle [rad]
    thrust.beta_rth = beta_rth_r2;
    
    % Azimuth angle [rad]
    thrust.alpha_rth = alpha_rth_r2;
    
elseif flag_rth_r2 == 0 && flag_rth == 0
    thrust.epsilon_rth = 0;
    thrust.beta_rth = 0;
    thrust.alpha_rth = 0;
    thrust.thrust_rth = 0;
    thrust.rth_flag_eps_r2 = 0;
end



% =========================================================================
% Inertial Acceleration - check this 
% =========================================================================
if flag_SRP
    
    % SRP Acceleration [m/s^2] - modify to add the SC/SUn distance!
    a_inertial = constants.P_Sun * Cr * A_m_SRP ;
    
    % Adimensional acceleration [DU/TU^2]
    thrust.epsilon_inertial = a_inertial * 10^(-3) * (  (constants.TU*constants.sec_day)^2  / constants.DU);
    
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


% =========================================================================
% Thrust at perigee and apogee
% ========================================================================
if flag_peri_apo
    T_adim = T * 10^(-3) * (  (constants.TU*constants.sec_day)^2  / constants.DU);
    thrust.T_adim = T_adim;
    thrust.flag_peri_apo = 1;
else
    T_adim = 0;
    thrust.T_adim = T_adim;
    thrust.flag_peri_apo = 0;
end


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

input.geopotential.J2 = 0;
input.geopotential.J3 = 0;
input.geopotential.J4 = 0;
input.geopotential.J5 = 0;

input.flag_3rd_Sun = 0;
input.flag_3rd_Moon = 0;

input.drag.CD = 0;

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
    
    drag.CD = 0;
    
    % Settings for numerical propagation
    input_num_int.J2 = 0;
    input_num_int.drag.CD = 0;
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
    [T_Numeric, Equin_Numeric] = ode113(@(t,x)equations_motion3(t, x, ...
        input_num_int, thrust, control, drag, constants.g0),...
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
    figure
    plot((T_Analytic(:)-T_Analytic(1)) * constants.TU /365.25, h_perigee_Analytic * constants.DU,'b','LineWidth',2)
    hold on
    plot((T_Analytic(:)-T_Analytic(1)) * constants.TU/365.25, h_apogee_Analytic*constants.DU,'r','LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('Height [km]')
    legend('Perigee height','Apogee height')
    title('Analytic')
    
    
    % -------------------------------------------------------------------------
    % Keplerian elements 
    % -------------------------------------------------------------------------
    
    % Eccentricity and inclination
    e_Analytic = sqrt(Equin_Analytic(:,2).^2 + Equin_Analytic(:,3).^2);
    incl_Analytic = 2 * atan(sqrt(Equin_Analytic(:,4).^2 + Equin_Analytic(:,5).^2));
    
    
    
    figure('units','normalized','outerposition',[0 0 0.6 0.6])
    subplot(2,2,1)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/365.25, Keplerian_Analytic(1,:) * constants.DU,'LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('Semimajor axis [km]')
    subplot(2,2,2)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/365.25, Keplerian_Analytic(2,:),'LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('Eccentricity')
    subplot(2,3,4)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/365.25, Keplerian_Analytic(3,:)*180/pi,'LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('Inclination [deg]')
    subplot(2,3,6)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/365.25, Keplerian_Analytic(5,:)*180/pi,'LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('\omega [deg]')
    subplot(2,3,5)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/365.25, Keplerian_Analytic(4,:)*180/pi,'LineWidth',2)
    grid on
    xlabel('Time [years]')
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
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/365.25, Equin_Analytic(:,1) * constants.DU,'LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('Semimajor axis [km]')
    subplot(2,2,2)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/365.25, Equin_Analytic(:,2),'LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('P_1')
    subplot(2,3,4)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/365.25, Equin_Analytic(:,3),'LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('P_2')
    subplot(2,3,6)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/365.25, Equin_Analytic(:,4)*180/pi,'LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('Q_1 [deg]')
    subplot(2,3,5)
    plot((T_Analytic(:)-T_Analytic(1))*constants.TU/365.25, Equin_Analytic(:,5)*180/pi,'LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('Q_2')
    
    
    
    
    % -------------------------------------------------------------------------
    % 3D plot
    % -------------------------------------------------------------------------
    figure
    % Initial trajectory
    plot_trajectory(Cartesian_Analytic(1:3,1)*constants.DU,...
        Cartesian_Analytic(4:6,1)*constants.DU/constants.TU,...
        2*pi*sqrt((Keplerian_Analytic(1,1))^3/constants.mu)*constants.TU,...
        constants.mu*constants.DU^3/constants.TU^2,'k',4);
    hold on
    % Final trajectory
    plot_trajectory(Cartesian_Analytic(1:3,end)*constants.DU,...
        Cartesian_Analytic(4:6,end)*constants.DU/constants.TU,...
        2*pi*sqrt((Keplerian_Analytic(1,end))^3/constants.mu)*constants.TU,...
        constants.mu*constants.DU^3/constants.TU^2,'b',4);
    hold on
    grid on
    quiver3(0,0,0,5*constants.DU,0,0,'Color','k','LineWidth',2)
    quiver3(0,0,0,0,5*constants.DU,0,'Color','k')
    quiver3(0,0,0,0,0,5*constants.DU,'Color','k')
      
    for i = 2 : 400 : length(Keplerian_Analytic) -1
        
        plot_trajectory(Cartesian_Analytic(1:3,i)*constants.DU,...
            Cartesian_Analytic(4:6,i)*constants.DU/constants.TU, ...
            2*pi*sqrt((Keplerian_Analytic(1,i))^3/constants.mu)*constants.TU,...
            constants.mu*constants.DU^3/constants.TU^2,[0.5 0.5 0.5],1);
    end
    hold on
    [x_s,y_s,z_s] = sphere();
    C = ones(size(z_s))*0;
    hSurface = surf( x_s, y_s,z_s ,C) ;
    set(hSurface,'FaceColor',[0 0 0],'FaceAlpha',0);
    grid on
    xlabel('x [km]','FontSize',12)
    ylabel('y [km]','FontSize',12)
    zlabel('z [km]','FontSize',12)
    legend('Initial','Final')
    axhand=gca;
    set(axhand,'FontName','Helvetica','FontSize',12)
    view(3)
    axis equal
    
    
    
    
    % -------------------------------------------------------------------------
    % 3D plot - AU
    % -------------------------------------------------------------------------
    figure
    % Initial trajectory
    plot_trajectory(Cartesian_Analytic(1:3,1),...
        Cartesian_Analytic(4:6,1),...
        2*pi*sqrt((Keplerian_Analytic(1,1))^3/constants.mu),...
        constants.mu,'k',4);
    hold on
    % Final trajectory
    plot_trajectory(Cartesian_Analytic(1:3,end),...
        Cartesian_Analytic(4:6,end),...
        2*pi*sqrt((Keplerian_Analytic(1,end))^3/constants.mu),...
        constants.mu,'b',4);
    hold on
    grid on
    quiver3(0,0,0,5,0,0,'Color','k','LineWidth',2)
    quiver3(0,0,0,0,5,0,'Color','k')
    quiver3(0,0,0,0,0,5,'Color','k')
      
    for i = 2 : 400 : length(Keplerian_Analytic) -1
        
        plot_trajectory(Cartesian_Analytic(1:3,i),...
            Cartesian_Analytic(4:6,i), ...
            2*pi*sqrt((Keplerian_Analytic(1,i))^3/constants.mu),...
            constants.mu,[0.5 0.5 0.5],1);
    end
    hold on
    [x_s,y_s,z_s] = sphere();
    C = ones(size(z_s))*0;
%     hSurface = surf( x_s, y_s,z_s ,C) ;
    set(hSurface,'FaceColor',[0 0 0],'FaceAlpha',0);
    grid on
    xlabel('x [AU]','FontSize',12)
    ylabel('y [AU]','FontSize',12)
    zlabel('z [AU]','FontSize',12)
    legend('Initial','Final')
    axhand=gca;
    set(axhand,'FontName','Helvetica','FontSize',12)
    view(3)
    axis equal

    
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
    plot((T_Numeric(:)-T_Numeric(1)) * constants.TU /365.25, h_perigee_Numeric * constants.DU,'r')
    hold on
    plot((T_Numeric(:)-T_Numeric(1))*constants.TU/365.25,   h_apogee_Numeric*constants.DU,'r')
    grid on
    xlabel('Time [years]')
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
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU / 365.25, h_perigee_Numeric * constants.DU,'r')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU / 365.25, h_perigee_Analytic * constants.DU,'b','LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('Perigee height [km]')
    legend('Numerical','Averaged analytical')
    
    subplot(1,2,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/365.25,   h_apogee_Numeric*constants.DU,'r')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/365.25, h_apogee_Analytic*constants.DU,'b','LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('Apogee height [km]')
    legend('Numerical','Averaged analytical')
    
    
    
    % -------------------------------------------------------------------------
    % Equinoctial elements for numerical and averaged analytical solution: a,
    % P1, P2 and Q1
    % -------------------------------------------------------------------------
    figure
    title('Equinoctial elements')
    subplot(2,2,1)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/365.25, Equin_Analytic(:,1)*constants.DU,'b')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU / 365.25,Equin_Numeric(1:steps_num:end,1)*constants.DU,'r')
    xlabel('Time [years]')
    ylabel('a [km]')
    legend('Averaged analytical','Numerical')
    
    subplot(2,2,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU / 365.25, Equin_Analytic(:,2),'b')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/365.25, Equin_Numeric(1:steps_num:end,2),'r')
    xlabel('Time [years]')
    ylabel('P_1')
    legend('Averaged analytical','Numerical')
    
    subplot(2,2,3)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, Equin_Analytic(:,3),'b')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, Equin_Numeric(1:steps_num:end,3),'r')
    xlabel('Time [years]')
    ylabel('P_2')
    legend('Averaged analytical','Numerical')
    
    subplot(2,2,4)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, Equin_Analytic(:,4),'b')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, Equin_Numeric(1:steps_num:end,4),'r')
    xlabel('Time [years]')
    ylabel('Q_1')
    legend('Averaged analytical','Numerical')
    
    
    figure
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, Equin_Analytic(:,5),'b')
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, Equin_Numeric(1:steps_num:end,5),'r')
    xlabel('Time [years]')
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
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU /365.25,Equin_Numeric(1:steps_num:end,1)*constants.DU,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/365.25, Equin_Analytic(:,1)*constants.DU,'b','LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('a [km]')
    legend('Numerical','Averaged analytical')
    
    subplot(2,1,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/365.25, e_Numeric,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU /365.25, e_Analytic,'b','LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('e')
    legend('Numerical','Averaged analytical')
    
    % -------------------------------------------------------------------------
    % Inclination, RAAN and perigee argument
    figure
    subplot(2,2,1)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, incl_Numeric*180/pi,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, incl_Analytic*180/pi,'b','LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('i [deg]')
    legend('Numerical','Averaged analytical')
    
    subplot(2,2,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, Keplerian_Numeric(4,:)*180/pi,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, Keplerian_Analytic(4,:)*180/pi,'b','LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('\Omega [deg]')
    legend('Numerical','Averaged analytical')
    
    subplot(2,2,3)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, Keplerian_Numeric(5,:)*180/pi,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, Keplerian_Analytic(5,:)*180/pi,'b','LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('\omega [deg]')
    legend('Numerical','Averaged analytical')
    
    
    
        figure('units','normalized','outerposition',[0 0 0.6 0.6])
    subplot(2,2,1)
       plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU /365.25,Equin_Numeric(1:steps_num:end,1)*constants.DU,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/365.25, Equin_Analytic(:,1)*constants.DU,'b','LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('a [km]')
    legend('Numerical Propagation','Averaged Analytical Propagation')
    axhand=gca;
set(axhand,'FontSize',14)
    subplot(2,2,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/365.25, e_Numeric,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU /365.25, e_Analytic,'b','LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('e')
        axhand=gca;
set(axhand,'FontSize',14)
%     legend('Numerical','Averaged analytical')
    subplot(2,3,4)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, incl_Numeric*180/pi,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, incl_Analytic*180/pi,'b','LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('i [deg]')
        axhand=gca;
set(axhand,'FontSize',14)
%     legend('Numerical','Averaged analytical')
    subplot(2,3,6)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, Keplerian_Numeric(4,:)*180/pi,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, Keplerian_Analytic(4,:)*180/pi,'b','LineWidth',2)
    grid on
    xlabel('Time [years]')
    ylabel('\Omega [deg]')
        axhand=gca;
set(axhand,'FontSize',14)
%     legend('Numerical','Averaged analytical')
    subplot(2,3,5)
 plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, (Keplerian_Analytic(5,:)-Keplerian_Numeric(5,:))*180/pi,'k','LineWidth',2)
plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, Keplerian_Numeric(5,:)*180/pi,'r','LineWidth',2)
    hold on
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/constants.sec_day, Keplerian_Analytic(5,:)*180/pi,'b','LineWidth',2)
    grid on
    xlabel('Time [years]')
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
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/365.25, (Equin_Analytic(:,1)-Equin_Numeric(1:steps_num:end,1))*constants.DU,'k','LineWidth',2)
    xlabel('Time [years]')
    ylabel('\Delta a [km]')
    subplot(2,1,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU /365.25, e_Analytic-e_Numeric,'k','LineWidth',2)
    xlabel('Time [years]')
    ylabel('\Delta e')
    
    
    figure
    subplot(2,2,1)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, (incl_Analytic-incl_Numeric)*180/pi,'k','LineWidth',2)
   xlabel('Time [years]')
    ylabel('\Delta i [deg]')
    subplot(2,2,2)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, (Keplerian_Analytic(4,:)-Keplerian_Numeric(4,:))*180/pi,'k','LineWidth',2)
   xlabel('Time [years]')
    ylabel('\Delta \Omega [deg]')
    subplot(2,2,3)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, (Keplerian_Analytic(5,:)-Keplerian_Numeric(5,:))*180/pi,'k','LineWidth',2)
   xlabel('Time [years]')
    ylabel('\Delta \omega [deg]')
    
    
    
    
    
    figure('units','normalized','outerposition',[0 0 0.6 0.6])
    subplot(2,2,1)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU/365.25, (Equin_Analytic(:,1)-Equin_Numeric(1:steps_num:end,1))*constants.DU,'k','LineWidth',2)
    xlabel('Time [years]')
    ylabel('\Delta a [km]')
    subplot(2,2,2)
   plot((T_Numeric(1:steps_num:end)-T_Numeric(1))*constants.TU /365.25, e_Analytic-e_Numeric,'k','LineWidth',2)
   xlabel('Time [years]')
    ylabel('\Delta e')
    subplot(2,3,4)
   plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, (incl_Analytic-incl_Numeric)*180/pi,'k','LineWidth',2)
  xlabel('Time [years]')
    ylabel('\Delta i [deg]')
    subplot(2,3,6)
    plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/365.25, (Keplerian_Analytic(4,:)-Keplerian_Numeric(4,:))*180/pi,'k','LineWidth',2)
   xlabel('Time [years]')
    ylabel('\Delta \Omega [deg]')
    subplot(2,3,5)
 plot((T_Numeric(1:steps_num:end)-T_Numeric(1)) * constants.TU/364.25, (Keplerian_Analytic(5,:)-Keplerian_Numeric(5,:))*180/pi,'k','LineWidth',2)
    xlabel('Time [years]')
    ylabel('\Delta \omega [deg]')
    
    
    
    
    
    
end



