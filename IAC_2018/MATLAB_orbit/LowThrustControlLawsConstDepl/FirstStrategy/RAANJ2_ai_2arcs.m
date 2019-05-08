function [DeltaV, x_alpha, beta, ToF_1st, ToF_2nd] = RAANJ2_ai_2arcs(initial_orbit, final_orbit, ToF_total, ...
    f, constants, plot_flag)




% Input: initial_orbit -> structure with the initial orbital elements
%        final_orbit   -> structure with the final orbital elements
%        ToF_total     -> total time of fligth
%        f             -> low-thrust acceleration
%        constants     -> structure with constants (J2, gravitational par.)
%        plot_flag     ->

% Output: DeltaV -> cost of the transfer
%         x_alpha -> semiamplitude of the thrust arcs for the second phase
%         beta    -> elevation angle during second phase
%         ToF1   -> time of flight of the first phase
%         ToF2   -> time of flight of the second phase

% Marilena Di Carlo, 2016
% marilena.di-carlo@strath.ac.uk

%%

a_initial    = initial_orbit.a;
i_initial    = initial_orbit.incl;
RAAN_initial = initial_orbit.RAAN;

a_final    = final_orbit.a;
i_final    = final_orbit.incl;
RAAN_final = final_orbit.RAAN;

R_Earth = constants.R_Earth;
J2      = constants.J2;
mu      = constants.mu;

%% Compute alpha to realise transfer in given time of flight

% Variation of RAAN on the initial orbit
RAAN_dot = (1.5 * sqrt(mu) * J2 * R_Earth^2 * cos(i_initial) *  a_initial^(-7/2));

% Initialise n to zero
n = 0;

% Equations are singular when i_final is approximately equal to i_initial.
% The first if condition consider the case in which the initial and final
% inclination have a difference greater than 0.5 deg
if abs(i_final - i_initial) > 0.5*pi/180  && a_final ~= a_initial
    
    % Parameter for ToF_1st
    p1 = ((RAAN_final - RAAN_initial) ) / RAAN_dot;
    
    % Here a for cycle is required to increase every time the value of "n"
    % that appears in ff
    for ll =  1 : 10
        
        
        k2 = - 4 * log(a_final/a_initial) / (i_final - i_initial);
        
        k3 = exp(k2 * i_final) * (k2 * cos(i_final) + sin(i_final)) - ...
            exp(k2 * i_initial) * (k2 * cos(i_initial) + sin(i_initial));
        
        k1_ = -3/4 * mu * pi * J2 * R_Earth^2 / (a_initial^4 * f) * ...
            exp(4 * log(a_final/a_initial) * i_initial / (i_final - i_initial));
        
        p2 = k1_ * k3 / (1+k2^2) / RAAN_dot;
        
        p3 = pi * sqrt(mu) / (4*f) * (2/sqrt(a_initial) - 2/sqrt(a_final)) * ...
            1/ log(a_final/a_initial);
        
        p2_ = p2 / (2 * (i_final - i_initial));
        
        
        % Equation to solve:
        ff = @(xa)( -p1+ p2_ * sqrt((log(a_final/a_initial))^2 * sin(xa)^2 + ...
            4 * (i_final - i_initial)^2 * xa^2) / ...
            (xa * sin(xa))  +  2*(n)*pi/RAAN_dot+ ...
            p3 * sqrt((log(a_final/a_initial))^2 * sin(xa)^2 + ...
            4 * (i_final - i_initial)^2 * xa^2) / ...
            (xa * sin(xa)) -...
            ToF_total);
        
        
        % Solve equation
        [x_alpha, ~, exitflag] = fzero(ff,0.5);
        
        % Time of flight second phase
        ToF_2nd_foreseen = p3 * sqrt((log(a_final/a_initial))^2 * sin(x_alpha)^2 + ...
            4 * (i_final - i_initial)^2 * x_alpha^2) / ...
            (x_alpha * sin(x_alpha));
        
        % Time of flight first phase
        ToF_1st =   -p1 + p2_ * sqrt((log(a_final/a_initial))^2 * sin(x_alpha)^2 + ...
            4 * (i_final - i_initial)^2 * x_alpha^2) / ...
            (x_alpha * sin(x_alpha)) + 2 * (n) * pi / RAAN_dot;
        
        
        % Increase n
        n = n + 1;
        
        % Break cycle and accept solution when:
        % Time of flight of first phase is greater than zero
        % Time of flight of first phase is such that there is not a
        % variation of RAAN in the first phase that is greater than 2*pi
        % The semi-amplitude of the thrust arc (alpha) is greater than 0 but
        % lower than 90 deg
        if ToF_1st > 0  && abs(ToF_1st)*RAAN_dot < 4*pi && x_alpha>0 && x_alpha<pi/2
            break
        end
        
    end
    
    % If for more than 10 cycle the above process did not work, then exit
    % with NaN results
    if ll == 10 && (ToF_1st<0 || ToF_2nd_foreseen<0)
        DeltaV = NaN;
        x_alpha = NaN;
        beta = NaN;
        ToF_2nd = NaN;
        return
    end
    

 % If initial and final inclination are almost the same   
elseif abs(i_final - i_initial) < 0.5*pi/180 && a_initial ~= a_final
    
    n = 0;
    exitflag = 0;
    for ll = 1 : 10

        
        k4 = 3/32 * pi * mu * J2 * R_Earth^2 * cos(i_initial) / f * ...
            (1/(a_final^4) - 1/(a_initial^4)) ;
        
        x_alpha = (pi/(2 * f) * (sqrt(mu/a_initial) - sqrt(mu/a_final)) + k4/RAAN_dot) / ...
            (ToF_total + (RAAN_final - RAAN_initial) / RAAN_dot - 2 * (n) * pi / RAAN_dot);
            
        ToF_2nd_foreseen = pi / (2 * f * x_alpha) * (sqrt(mu/a_initial) - sqrt(mu/a_final));
        
        ToF_1st = -(RAAN_final - k4/x_alpha - RAAN_initial)/ RAAN_dot + 2 * (n) * pi / RAAN_dot;

        
        n = n + 1;
        if ToF_1st > 0 && abs(ToF_1st)*RAAN_dot < 4*pi && abs(x_alpha) < pi/2
            exitflag = 1;
            break
        end

    end

  % If there is no variation of semimajor axis  
elseif abs(i_final - i_initial)>0.5*pi/180 && a_initial == a_final
    
    n = 0;
    exitflag = 0;
    for ll = 1 : 100
        
        num = 3/4 * pi * mu * J2 * R_Earth^2 * (sin(i_final) - sin(i_initial)) / ...
            (a_initial^4 * f * RAAN_dot) + ...
            pi/(2 * f) * sqrt(mu/a_initial) * (i_final - i_initial);
        
        den = ToF_total - (RAAN_final - RAAN_initial) / RAAN_dot +  2 * (n) * pi / RAAN_dot;
        
        x_alpha = asin(num/den);
        
        k1 = -3/4 * pi * mu * J2 * R_Earth^2 / (a_initial^4 * f * sin(x_alpha));
        
        ToF_1st = ( RAAN_final - k1 * (sin(i_final) - sin(i_initial)) - RAAN_initial) / RAAN_dot + ...
            - 2 * (n) * pi / RAAN_dot ;
        
        ToF_2nd_foreseen = pi/(2*f) * sqrt(mu/a_initial) * (i_final - i_initial) / sin(x_alpha);
        
        n = n + 1;
        if ToF_1st > 0 && abs(ToF_1st)*RAAN_dot < 2*pi
            
            
            exitflag = 1;
            break
        end
        
    end
    
    
end

% If there was a solution to the problem
if exitflag == 1
    
    % Delta_V [km/s] and ToF [h] required for the transfer of the second
    % phase
    [DeltaV, ToF_2nd, DeltaRAAN_2nd, beta] = a_i_circular_2arcs(initial_orbit, final_orbit, ...
        x_alpha, f,...
        constants, plot_flag);
    
else
    DeltaV = NaN;
    return
end


if abs(ToF_2nd - ToF_2nd_foreseen)>1e-6
    keyboard
end


%% 1st phase

% If the sum of the time of flight gives an error
if abs(ToF_1st + ToF_2nd - ToF_total) > 1e-6
    warning('Wrong solution: sum of ToF does not correspond to total ToF')
    keyboard
end


if plot_flag
    
    t = linspace(0, ToF_1st, 1000);
    
    for i = 1 : length(t)
        RAAN_1st(i) = RAAN_initial -1.5 * sqrt(mu) * J2 * R_Earth^2 * ...
            cos(i_initial) * a_initial^(-7/2) * t(i);
    end
    
    figure
    plot(t/86400, mod(RAAN_1st, 2*pi)*180/pi, 'LineWidth',2)
    grid on
    xlabel('Time [days]')
    ylabel('RAAN [deg]')
    
    %     keyboard
    
    [DeltaV, ToF_2nd, DeltaRAAN_2nd, beta] = a_i_circular_2arcs(initial_orbit, final_orbit, ...
        x_alpha, f,...
        constants, plot_flag);
    
    %     keyboard
end

