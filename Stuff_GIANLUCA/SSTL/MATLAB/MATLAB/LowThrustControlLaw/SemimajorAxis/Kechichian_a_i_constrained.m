function [Delta_V, ToF_hours] = Kechichian_a_i_constrained(a_initial, a_final, a_limit, i_initial, i_final, f, mu, plot_flag)

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

% Marilena Di Carlo
% marilena.di-carlo@strath.ac.uk


% Velocity on initial orbit [km/s]
v_initial = sqrt(mu / a_initial);

% Velocity on final orbit [km/s]
v_final = sqrt(mu / a_final);

% Limit velocity [km/s]
v_limit = sqrt(mu / a_limit);

% Inclination difference
Delta_i0 = (i_final - i_initial);

% =========================================================================
% Initial yaw angle

% Initial yaw angle
sin_beta0 = (v_limit / v_initial);
cos_beta0 = sqrt(1 - sin_beta0^2);
beta0 = atan2(sin_beta0, cos_beta0);
% =========================================================================
% keyboard
% Constraint entry time
t1 = v_initial * cos_beta0 / f;

% Constraint exit time
num = sqrt((v_final - v_limit) * (v_final + v_limit));
den = v_limit;
tan_term = atan2(num, den);
t2 = pi * v_limit / (2*f) * (i_final - 2/pi * (pi/2 - beta0) - 2/pi * tan_term) + (v_initial * cos_beta0 / f);

% DeltaV required for the transfer [km/s]
% Kechichian 2010, Equation 6.20


% Time of flight [s]
ToF = t2 + num / f;

% Time of flight [h]
ToF_hours = ToF / 86400;

% Time vector
t = linspace(0, ToF, 1000);

incl_t1 = i_initial + 2/pi * (atan2((f * t1 - v_initial * cos(beta0)), v_initial * sin(beta0)) + pi/2 - beta0);
incl_t2 = incl_t1 + 2 * f / (pi * v_limit) * (t2 - t1);
% Yaw angle at time t - sin


% Equation 6.17 Kechichian 2010


% Initialization
cos_beta_t       = zeros(1,numel(t));
beta_t           = zeros(1,numel(t));
v_t              = v_limit * ones(1,numel(t));
Delta_i          = zeros(1, numel(t));
sin_Delta_i_term = zeros(1,numel(t));

% Compute time varying quantities
for index = 1 : numel(t)
    
    if t(index) < t1
        % Velocity in the interval [0, t1]
        v_t(index) = sqrt(v_initial^2 + f^2 * t(index)^2 - 2 * f * t(index) * v_initial * cos_beta0);
        
        % Yaw angle in the interval [0, t1]
        beta_t(index) = atan2(v_limit, v_initial * cos_beta0 - f * t(index));
        
        incl(index) = i_initial + 2/pi * (atan2((f * t(index) - v_initial * cos(beta0)), v_initial * sin(beta0)) + pi/2 - beta0);
        
    elseif t(index) > t1 && t(index) < t2
        beta_t(index) = beta_t(index-1);
        
        incl(index) = incl_t1 + 2 * f / (pi * v_limit) * (t(index) - t1);
        
    elseif t(index) > t2
        % Velocity in the interval [t2, TOF]
        v_t(index) = sqrt(v_limit^2 + f^2 * (t(index) - t2)^2 );
        
        % Yaw angle in the interval [t2, TOF]
        beta_t(index) = atan2(-v_limit, f * (t(index) - t2));
        
        incl(index) = i_initial +  2/pi * (pi/2 - beta0) + 2 * f /(pi * v_limit) * (t2 - v_initial * cos_beta0 / f) + ...
            2/pi * atan2(f * (t(index) - t2), v_limit);
        
    end
    
    % Inclination  at time t [rad]

end


              
% Velocity at t1
v_t1 = sqrt(v_initial^2 + f^2 * t1^2 - 2 * f * t1 * cos(beta0));

% =========================================================================
% DeltaV
% =========================================================================
% From t0 to t1
% Equation 43, paper:
DeltaV_1 = v_initial * cos(beta0) - v_initial * sin(beta0) / (tan((pi/2)*(incl_t1 - i_initial) + beta0));

% From t2 to tf - Equation 43
% DeltaV_2 = v_limit * cos(pi/2) + sqrt(v_final^2 - v_limit^2 * sin(pi/2)^2);
DeltaV_2 = v_limit * cos(pi/2) - v_limit * sin(pi/2) / (tan((pi/2)*(incl(end) - incl_t2) + pi/2));

% From t1 to t2
% DeltaV_12 = sqrt(v_limit^2 + 2 * v_limit^2 * cos(pi/2 * (incl_t2 - incl_t1)) + v_limit^2);
% v_limit * cos(pi/2)  - v_limit * sin(pi/2) / (tan((pi/2)*(incl_t2 - incl_t1) + pi/2))
% DeltaV2 = t2 *f + v_limit * cos(pi/2) + sqrt(v_limit^2 * cos(pi/2)^2 + v_final^2 - 4 * v_limit^2);
DeltaV_12 = pi * v_limit * (incl_t2 - incl_t1)/ 2 ;

Delta_V = DeltaV_1  + DeltaV_2 + DeltaV_12;

% keyboard

% Semimajor axis at time t [km]
a = mu ./ v_t.^2;



if plot_flag == 1
    % Plot results
    figure
    subplot(2,2,1)
    plot(t/86400, a,'LineWidth',2)
    hold on
    line([0 ToF_hours],[a_final a_final],'Color','r','LineWidth',2)
    xlabel('Time [days]')
    ylabel('Semimajor axis [km]')
    grid on
    subplot(2,2,2)
    plot(t/86400, incl*180/pi,'LineWidth',2)
    hold on
    line([0 ToF_hours],[i_final*180/pi i_final*180/pi],'Color','r','LineWidth',2)
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





