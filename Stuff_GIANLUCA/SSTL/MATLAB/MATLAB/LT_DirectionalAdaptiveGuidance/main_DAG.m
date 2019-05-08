% -------------------------------------------------------------------------
% Original Directional Adaptive Guidance (DAG)
% Reference: Ruggiero, Low-thrust manuevers for the efficient correction of
% orbital elements
% -------------------------------------------------------------------------
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

clear
close all


%% Constants - Earth transfer

% Gravitational parameter [km^3/s^2]
constants.mu_dim = 398600;
% Adimensional gravitational parameter
constants.mu = 1;
% Distance unit (Earth radius) [km]
constants.DU = 6378.136;
% Time unit [s]
constants.TU = 806.78;
% Adimensional g0
constants.g0 = 9.8 * 1e-3 * constants.TU^2 / constants.DU;
constants.g0_dim = 9.8;
% 
constants.sec_day = 86400;

%% Input

% Thrust [N]
engine.T_dim = 0.089;
% Specific impulse [s]
engine.Isp_dim = 1650;

% Spacecraft mass [kg]
m0 = 367;

% Initial orbital elements (semimajor axis in DU, angles in radians)
x0 = [24396/constants.DU ...
     0.7283 ...
     7*pi/180 ...
     0 ...
     0 ...
     0];

% Final targeted orbital elements
xf = [42164/constants.DU ...
      0 ...
      0*pi/180 ...
      0 ...
      0 ...
      0];

% Targeted orbital elements? 1 to target, 0 to leave unconstrained
target_flag = [1 1 1 0 0];

% Tolerance on manuever efficiency
tol = 0.6;

% Maximum time of flight [days]
ToF = 200;



%% 

% Adimensional engine parameters
engine.T   = engine.T_dim * 1e-3 / constants.DU * (constants.TU)^2;
engine.Isp = engine.Isp_dim / constants.TU;

% DAG parameters
DAG_parameters.x0          =  x0;
DAG_parameters.x_target    = xf;
DAG_parameters.target_flag = target_flag;
DAG_parameters.tol         = tol;

% Set options for numerical integration
options = odeset('RelTol',1e-8,'AbsTol',1e-10,...
        'events',@(t,x)events_stop_ode(t,x,DAG_parameters));

[t,x_final] = ode113(@(t,x)Gauss_eqns(t, x, DAG_parameters, engine, constants),...
    [0 ToF*constants.sec_day/constants.TU],[DAG_parameters.x0 m0],options);

% DeltaV [km/s]
DeltaV = engine.Isp_dim * constants.g0_dim * log(x_final(1,7)/x_final(end,7)) / 1000;

% Control history - find different way to compute this
f_LT_history = zeros(3,length(t));
for index = 1 : length(x_final)
   
    f_LT_history(:,index) = DAG_LT_acc(x_final(index,1:6), DAG_parameters, constants);
    
end

% Acceleration magnitude
f_LT.magnitude = (engine.T ./ x_final(:,7))' .* ...
    sqrt(f_LT_history(1,:).^2 + f_LT_history(2,:).^2 +  f_LT_history(3,:).^2);
% f_LT.magnitude(index) = norm(f_LT_history(:,index));

% Azimuth angle
f_LT.alpha = atan2(f_LT_history(1,:), f_LT_history(2,:));

% Elevation angle
f_LT.beta = atan2(f_LT_history(3,:), ...
    sqrt(f_LT_history(1,:).^2 + f_LT_history(2,:).^2));




%% Plot - orbital elements

figure
subplot(2,2,1)
plot(t*constants.TU/constants.sec_day, x_final(:,1) * constants.DU, 'LineWidth',2)
xlabel('Time [days]')
ylabel('a')
grid on
if target_flag(1)
    line([0 t(end)*constants.TU/constants.sec_day],[xf(1)*constants.DU xf(1)*constants.DU],...
        'Color','r','LineWidth',2)
end
subplot(2,2,2)
plot(t*constants.TU/constants.sec_day, x_final(:,2), 'LineWidth',2)
xlabel('Time [days]')
ylabel('e')
grid on
if target_flag(2)
    line([0 t(end)*constants.TU/constants.sec_day],[xf(2) xf(2)],...
        'Color','r','LineWidth',2)
end
subplot(2,3,4)
plot(t*constants.TU/constants.sec_day, x_final(:,3) * 180/pi, 'LineWidth',2)
xlabel('Time [days]')
ylabel('i [deg]')
grid on
if target_flag(3)
    line([0 t(end)*constants.TU/constants.sec_day],[xf(3)*180/pi xf(3)*180/pi],...
        'Color','r','LineWidth',2)
end
subplot(2,3,5)
plot(t*constants.TU/86400, mod(x_final(:,4), 2*pi) * 180/pi, 'LineWidth',2)
xlabel('Time [days]')
ylabel('\Omega [deg]')
grid on
if target_flag(4)
    line([0 t(end)*constants.TU/constants.sec_day],[xf(4)*180/pi xf(4)*180/pi],...
        'Color','r','LineWidth',2)
end
subplot(2,3,6)
plot(t*constants.TU/constants.sec_day, mod(x_final(:,5), 2*pi) * 180/pi, 'LineWidth',2)
xlabel('Time [days]')
ylabel('\omega [deg]')
grid on
if target_flag(5)
    line([0 t(end)*constants.TU/constants.sec_day],[xf(5)*180/pi xf(5)*180/pi],...
        'Color','r','LineWidth',2)
end

%% Plot - control

figure
subplot(3,1,1)
plot(t * constants.TU / constants.sec_day, f_LT.magnitude)
grid on
xlabel('Time [days]')
ylabel('f [m/s^2]')
subplot(3,1,2)
plot(t * constants.TU / constants.sec_day, (f_LT.alpha)*180/pi)
grid on
xlabel('Time [days]')
ylabel('\alpha [deg]')
subplot(3,1,3)
plot(t * constants.TU / constants.sec_day, (f_LT.beta)*180/pi)
grid on
xlabel('Time [days]')
ylabel('\beta [deg]')