function [Delta_V, ToF_hours, a] = Edelbaum_a_i(a_initial, a_final, i_initial, i_final, f, mu, plot_flag)

% Input: a_initial -> initial semimajor axis
%        a_final   -> final semimajor axis
%        i_initial -> initial inclinaiton
%        i_final   -> final inclination
%        f         -> acceleration [km/s^"]
%        mu        -> gravitational acceleration [km^3 /s^2]
%        plot_flag -> 1 to realize the plot, 0 otherwise

% Output: Delta_V   -> Delta V [km/s] to realize the transfer
%         ToF_hours -> time of flight in hours

% References: Kechichian, "Reformulation of Edelbaum's Low-Thrust Transfer
% Problem Using Optimal Control Theory", 1997;
% Kechichian, "Analytical Representations of Optimal Low-Thurst Transfer in
% Circular Orbit, 2010"

% Marilena Di Carlo, 2015
% marilena.di-carlo@strath.ac.uk



% Velocity on initial orbit [km/s]
v_initial = sqrt(mu / a_initial);

% Velocity on final orbit [km/s]
v_final = sqrt(mu / a_final);

% Inclination difference
Delta_i = (i_final - i_initial);



% =========================================================================
% Initial yaw angle

% Initial yaw angle - sin
sin_beta0 = sin(pi/2 * Delta_i);
% Initial yaw angle - cos
cos_beta0 = v_initial / v_final - cos(pi/2 * Delta_i);
% Initial yaw angle
beta0 = atan2(sin_beta0, cos_beta0);
% =========================================================================



% DeltaV required for the transfer [km/s]
% Kechichian 2010, Equation 6.20
if beta0 ~= 0 && beta0 ~= pi
    Delta_V = v_initial * cos(beta0) - v_initial * sin(beta0) / (tan(pi/2*Delta_i + beta0));
else
    Delta_V = sqrt(v_initial^2 - 2 * v_initial * v_final + v_final^2);
end


% Time of flight [s]
ToF = Delta_V / f;

% Time of flight [h]
ToF_hours = ToF / 86400;

% Time vector
t = linspace(0, ToF, 100);

% Yaw angle at time t - sin
sin_beta_t       = v_initial * sin(beta0);

% Equation 6.17 Kechichian 2010
cos_Delta_i_term = v_initial * sin(beta0);

% Initialization
cos_beta_t       = zeros(1,numel(t));
beta_t           = zeros(1,numel(t));
v_t              = zeros(1,numel(t));
Delta_i          = zeros(1, numel(t));
sin_Delta_i_term = zeros(1,numel(t));

% Compute time varying quantities
for index = 1 : numel(t)
    
    % Yaw angle at time t - cosine
    cos_beta_t(index) = v_initial * cos(beta0) - f  * t(index);
    
    % Yaw angle at time t [rad]
    beta_t(index) = atan2(sin_beta_t,cos_beta_t(index));
    
    % velocity at time t [km/s]
    v_t(index) = sqrt(v_initial^2 + f^2 * t(index)^2 - 2 * f * t(index) * v_initial * cos(beta0));
    
    % 
    sin_Delta_i_term(index) = f * t(index) - v_initial * cos(beta0);
    
    % Inclination difference at time t (wrt initial inclination) [rad]
    Delta_i(index) = 2/pi * (beta_t(index) - beta0);
end

% Inclination at time t [rad]
incl = Delta_i + i_initial;

% Semimajor axis at time t [km]
a = mu ./ v_t.^2;


if plot_flag == 1
    % Plot results
    figure
    subplot(2,2,1)
    plot(t/86400, a,'LineWidth',2)
    hold on
    line([0 ToF_hours],[a_final a_final],'Color','r')
    xlabel('Time [days]')
    ylabel('Semimajor axis [km]')
    grid on
    subplot(2,2,2)
    plot(t/86400, incl*180/pi,'LineWidth',2)
    hold on
    line([0 ToF_hours],[i_final*180/pi i_final*180/pi],'Color','r')
    xlabel('Time [days]')
    ylabel('Inclination [deg]')
    grid on
    subplot(2,2,4)
    plot(t/86400, abs(beta_t)*180/pi,'LineWidth',2)
    xlabel('Time [days]')
    ylabel('\beta [deg]')
    grid on
    subplot(2,2,3)
    plot(t/86400, v_t,'LineWidth',2)
    xlabel('Time [days]')
    ylabel('Velocity [km/s]')
    grid on
end





