function [Delta_V, ToF_hours] = Kechichian_a_RAAN(a_initial, a_final, RAAN_initial, RAAN_final, incl, f, mu, plot_flag)

% Input: a_initial
%        a_final
%        i_initial
%        i_final
%        f

% Output: Delta_V   -> Delta V [km/s] to realize the transfer
%         ToF_hours -> time of flight in hours


% References: Kechichian, "Reformulation of Edelbaum's Low-Thrust Transfer
% Problem Using Optimal Control Theory", 1997;
% Kechichian, "Analytical Representations of Optimal Low-Thurst Transfer in
% Circular Orbit, 2010"

% Marilena Di Carlo, 2015



% Velocity on initial orbit [km/s]
v_initial = sqrt(mu / a_initial);

% Velocity on final orbit [km/s]
v_final = sqrt(mu / a_final);

% Right ascension difference [rad]
Delta_RAAN = (RAAN_final - RAAN_initial);
% keyboard
% =========================================================================
% Initial yaw angle

% Initial yaw angle - sin
sin_beta0 = sin(pi/2 * Delta_RAAN * sin(incl));
% Initial yaw angle - cos
cos_beta0 = v_initial / v_final - cos(pi/2 * Delta_RAAN * sin(incl));
% Initial yaw angle
beta0 = atan2(sin_beta0, cos_beta0);
% =========================================================================

% DeltaV required for the transfer [km/s]
Delta_V = sqrt( v_initial^2 - 2 * v_initial * v_final * cos(pi/2 * sin(incl) * Delta_RAAN) + v_final^2);

% Time of flight [s]
ToF = Delta_V / f;

% Time of flight [h]
ToF_hours = ToF / 86400;

% Time vector
t = linspace(0, ToF, 100);

% Yaw angle at time t - sin
sin_beta_t       = v_initial * sin(beta0);

% Equation 6.79 Kechichian 2010
cos_Delta_RAAN_term = sin(beta0);

% Initialization
cos_beta_t       = zeros(1,numel(t));
beta_t           = zeros(1,numel(t));
v_t              = zeros(1,numel(t));
Delta_RAAN          = zeros(1, numel(t));
sin_Delta_RAAN_term = zeros(1,numel(t));

% Compute time varying quantities
for index = 1 : numel(t)
    
    % Yaw angle at time t - cosine
    cos_beta_t(index) = v_initial * cos(beta0) - f  * t(index);
    
    % Yaw angle at time t [rad]
    beta_t(index) = atan2(sin_beta_t, cos_beta_t(index));
    
    % velocity at time t [km/s]
    v_t(index) = sqrt(v_initial^2 + f^2 * t(index)^2 - 2 * f * t(index) * v_initial * cos(beta0));
    
    % 
    sin_Delta_RAAN_term(index) = (f * t(index) - v_initial * cos(beta0))/v_initial;
    
    % Inclination difference at time t (wrt initial inclination) [rad]
%     Delta_RAAN(index) = 2/(pi * sin(incl)) * (pi/2 - beta0 + mod(atan2(sin_Delta_RAAN_term(index), cos_Delta_RAAN_term), 2*pi));
    Delta_RAAN(index) = 2/(pi * sin(incl)) * (pi/2 - beta0 + atan(sin_Delta_RAAN_term(index)/ cos_Delta_RAAN_term));

end

% RAAN at time t [rad]
Delta_RAAN = Delta_RAAN - Delta_RAAN(1);
RAAN = Delta_RAAN + RAAN_initial;

% Semimajor axis at time t [km]
a = mu ./ v_t.^2;


if plot_flag 

% Plot results
figure
subplot(2,2,1)
plot(t/86400, a,'LineWidth',2)
hold on
line([0 ToF_hours],[a_final a_final],'Color','r')
xlabel('Time [days]')
ylabel('Semimajor axis [km]')
% legend('','Target')
grid on
subplot(2,2,2)
plot(t/86400, RAAN*180/pi,'LineWidth',2)
hold on
line([0 ToF_hours],[RAAN_final*180/pi RAAN_final*180/pi],'Color','r')
xlabel('Time [days]')
ylabel('\Omega [deg]')
% legend('','Target')
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





