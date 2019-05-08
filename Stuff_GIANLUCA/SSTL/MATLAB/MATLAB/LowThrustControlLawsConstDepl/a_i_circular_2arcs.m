function [DeltaV, ToF, RAAN_ToF, beta] = a_i_circular_2arcs(initial_orbit, final_orbit, ...
    alpha, f, constants, plot_flag)

% Variation of semimajor axis and inclination using two thrust arc per
% revolution and constant elevation angle. 

% Input: intial_orbit -> structure with orbital elements of the initial
%                        orbit
%        final_orbit  -> structure with orbital elements of the final orbit
%        alpha        -> length thrust arcs
%        f            -> acceleration low thrust engine
%        constants    -> structure containing constants for the problem
%        plot_flag    ->

% Output: DeltaV   -> cost of the transfer
%         ToF      -> time of flight required for the transfer
%         RAAN_ToF -> final RAAN at the end of the transfer
%         beta     -> elevation angle to realise the transfer

% Reference: M. Di Carlo, A. Ricciardi, M. Vasile, "Optimised Constellation
% Deployment using Low-Thrust Propulsion", SPACE 2016

% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk
% 2016


%% Initialisation

% Initial orbital elements
a_initial    = initial_orbit.a;
i_initial    = initial_orbit.incl;
RAAN_0       = initial_orbit.RAAN;

% Final orbital elements
a_final    = final_orbit.a;
i_final    = final_orbit.incl;

% Constants
R_Earth = constants.R_Earth;
J2      = constants.J2;
mu      = constants.mu;

%% 

% Constant elevation angle:
beta = atan(2 * alpha * (i_final - i_initial) / sin(alpha) / log(a_final/a_initial));

% Time of flight:
% Special case with no variation of semimajor axis
if a_initial == a_final
    ToF = pi / (2 * f) * sqrt(mu/a_initial) * (i_final - i_initial) / (sin(alpha));
% Special case with no variation of inclination
elseif abs(i_final - i_initial) < 0.5 * pi/180 
    ToF = pi * sqrt(mu) / (4 * f  * alpha) * (2/sqrt(a_initial) - 2/sqrt(a_final));
% General case:
else
    ToF = pi * sqrt(mu) / (4 * f * cos(beta) * alpha) * (2/sqrt(a_initial) - 2/sqrt(a_final));
end

% DeltaV
DeltaV = sqrt(mu) / (2 * cos(beta))* (2/sqrt(a_initial) - 2/sqrt(a_final));

% Plot?
if plot_flag
    
    % Vector of time
    t = linspace(0, ToF, 1000);
    
    % Velocity and semimajor axis during the transfer
    V_t = sqrt(mu/a_initial) - 2 * f * cos(beta) * alpha * t / pi;
    a_t = mu./V_t.^2;
    
    % Variation of inclination during the transfer
    if a_initial == a_final
        i_t = i_initial + 2 * f / pi * sqrt(a_initial/mu) * sin(alpha) * t;
    else
        i_t = i_initial + tan(beta) * sin(alpha) / alpha * log(a_t/a_initial)/2;
    end
    
    figure
    subplot(1,2,1)
    plot(t/86400, a_t,'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('Semimajor axis [km]')
    subplot(1,2,2)
    plot(t/86400, i_t*180/pi,'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('Inclination [deg]')
    
    % General case: variation of semimajor axis AND inclination
    if beta ~= 0 && a_final ~= a_initial
        
        % Variation of right ascension
        n = 4;
        for ti = 1 : length(t)
            [xx, ww] = GaussQuad(n, 0, t(ti));
            for i = 1 : n
                V_t(i) = sqrt(mu/a_initial) - 2 * f * cos(beta) * alpha * xx(i) / pi;
                a_t(i) = mu./V_t(i).^2;
                i_t(i) = i_initial + tan(beta) * sin(alpha) / alpha * log(a_t(i)/a_initial)/2;         
                ff(i) = -3/2 * sqrt(mu) * J2 * R_Earth^2 * cos(i_t(i)) * a_t(i)^(-7/2);
            end
            RAAN_t(ti) = mod(RAAN_0  + dot(ff,ww), 2*pi);
        end
        
    % No variation of semimajor axis - time variation of right ascension can be computed analytically   
    elseif beta == pi/2

        for ti = 1 : length(t)
            RAAN_t(ti) = mod(RAAN_0 -3/4 * pi * mu * J2 * R_Earth^2 / ...
                (a_initial^4 * f * sin(alpha)) * ...
                (sin(i_t(ti)) - sin(i_initial)), 2*pi);
        end
    
    % No variation of inclination - time var of RAAN computed analytically   
    else
        
        for ti = 1 : length(t)
            RAAN_t(ti) = mod(RAAN_0 + ...
                3/32 * pi * mu/(f*alpha) * J2 * R_Earth^2 * cos(i_initial) * ...
                (1/a_t(ti)^4 - 1/a_initial^4), 2*pi);
        end
        
    end
    
    figure
    plot(t/86400, RAAN_t*180/pi,'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('RAAN [deg]')
    
end


% Final value of right ascension at the end of the transfer 
if beta ~= 0  && a_final ~= a_initial
    
    % ---------------------------------------------------------------------
    % Numerical integration
    % ---------------------------------------------------------------------
%     n = 4;
%     [xx, ww] = GaussQuad(n, 0, ToF);
%     for i = 1 : n
%         V_tof(i) = sqrt(mu/a_initial) - 2 * f * cos(beta) * alpha * xx(i) / pi;
%         a_tof(i) = mu./V_tof(i).^2;
%         i_tof(i) = i_initial + tan(beta) * sin(alpha) / alpha * log(a_tof(i)/a_initial)/2;       
%         ff_tof(i) = -3/2 * sqrt(mu) * J2 * R_Earth^2 * cos(i_tof(i)) * a_tof(i)^(-7/2);
%     end
%     RAAN_ToF = mod(RAAN_0 + dot(ff_tof,ww),2*pi);
    
    
    % ---------------------------------------------------------------------
    % Analytic integration
    % ---------------------------------------------------------------------
    k2 = - 4 * log(a_final/a_initial) / (i_final - i_initial);
    
    k3 = exp(k2 * i_final) .* (k2 * cos(i_final) + sin(i_final)) - ...
        exp(k2 * i_initial) * (k2 * cos(i_initial) + sin(i_initial));
    
    k1 = -3/4 * mu * pi *J2 * R_Earth^2 / ...
        (a_initial^4 * f * sin(beta) * sin(alpha)) * ...
        exp(4 * log(a_final/a_initial) * i_initial / (i_final - i_initial));
    
    % Variation of RAAN during the second phase
    RAAN_ToF = mod(RAAN_0 + k1 / (1+k2^2) * k3, 2*pi);
    

    
elseif beta == pi/2
    
    RAAN_ToF = mod(RAAN_0 -3/4 * pi * mu * J2 * R_Earth^2 / ...
                                (a_initial^4 * f * sin(alpha)) * ...
                                (sin(i_final) - sin(i_initial)), 2*pi);
                            
else

%     RAAN_ToF = mod(RAAN_0 + 3/16 * mu/f * J2 * R_Earth^2 * cos(i_initial) * ...
%                                   (1/a_final^4 - 1/a_initial^4), 2*pi);
                              
    RAAN_ToF = mod(RAAN_0 + 3/32 * pi * mu/(f*alpha) * J2 * R_Earth^2 * cos(i_initial) * ...
      (1/a_final^4 - 1/a_initial^4), 2*pi);
        
    
end






