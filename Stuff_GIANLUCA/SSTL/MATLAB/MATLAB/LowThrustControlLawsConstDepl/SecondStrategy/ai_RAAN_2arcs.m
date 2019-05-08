function [best_solution, all_solutions] = ai_RAAN_2arcs(initial_orbit, final_orbit, ...
    ToF_total, f_direction, f, constants)


% Input: initial_orbit -> structure containing the initial orbital elements
%        final_orbit   -> structure containing the final orbital elements
%        ToF_total     -> total time of flight available for the transfer
%        f_direction   -> direction of the out of plane thrust (+1/-1)
%        f             -> low-thrust acceleration
%        constants     -> structure of constants

% Output: best_solution -> structure containing data about the minimum
%                          deltaV solution
%         all_solutions -> structure containing data about all the
%                          solutions with different ToF

% Marilena Di Carlo, 2016
% marilena.di-carlo@strath.ac.uk

%%

% Initial orbital elements
a_initial    = initial_orbit.a;
i_initial    = initial_orbit.incl;
RAAN_initial = initial_orbit.RAAN;

% Final orbital elements
a_final    = final_orbit.a;
i_final    = final_orbit.incl;
RAAN_final = final_orbit.RAAN;

% Constants
R_Earth = constants.R_Earth;
J2      = constants.J2;
mu      = constants.mu;


%%

% Initial and final velocities
V_initial = sqrt(mu / a_initial);
V_final   = sqrt(mu / a_final);


% minimum edelbaum time to go from initial to final semimajor axis. This is
% used as lower limit for the time of the first phase 
min_ToF_ed = sqrt(V_final^2 + V_initial^2 - 2 * V_final * V_initial * cos(pi/2*(i_final - i_initial)) )  / f/86400;

% Possible values for the ToF of the first phase
T1st =linspace(min_ToF_ed, ToF_total/86400, 1000)*86400;

% T1st = 5*86400;
T2nd = ToF_total - T1st;

for index = 1 : length(T1st)
    
    
    %% 1st part: change of a and i in ToF1st
    
    % Compute alpha to achieve transfer in given ToF
    ff = @(xa)(pi / (2 * f * xa) * (V_initial - V_final) * ...
        sqrt(1 + (4 * xa^2) * (i_final - i_initial)^2 / (sin(xa) * log(a_final / a_initial))^2) - ...
        T1st(index));
    
    [alpha1(index), fval, exitflag] = fzero(ff,1.5);
    
    if alpha1(index) > pi/2
        exitflag = 0;
    end
    
    if exitflag == 1
        
        % Delta_V [km/s] and ToF [h] required for the transfer
        plot_flag = 0;
        [DeltaV1st(index), ToF, RAAN_ToF, betaa(index)] = a_i_circular_2arcs(initial_orbit, final_orbit, ...
            alpha1(index), f, constants, plot_flag);
        
    else
        alpha1(index) = NaN;
        alpha2(index) = NaN;
        DeltaV1st(index) = NaN;
        DeltaV2nd(index) = NaN;
        DeltaV_TOT(index) = NaN;
        betaa(index) = NaN;
        continue
    end
    
    
    %% 2nd phase: change right ascension to the final value
    
    
    % Variation of RAAN due to J2
    RAAN_dot = 3/2* sqrt(mu) * J2 * R_Earth^2 * cos(i_final) * ...
        a_final^(-7/2);
    
    % Max variation of RAAN that can be obtained in given ToF
    DeltaRAAN_max = (f_direction * 2 * f / (pi * sin(i_final)) * sqrt(a_final/mu)  * T2nd(index)+ ...
        - RAAN_dot * T2nd(index));
    
    
    n = 0;
    exitflag = 1;
    for kk = 1 : 2
        
        
        % Compute length thrust arc second phase
        num = (RAAN_final - RAAN_ToF + 1 *  2 * n * pi)/T2nd(index) + RAAN_dot ;
        den = f_direction * 2 * f / (pi * sin(i_final)) * sqrt(a_final/mu) ;
        
        alpha2(index) = asin(num/den);
             
        n = n+1;
        
        if isreal(alpha2(index)) && alpha2(index) >= 0
            break
        end
        
        %
% %         keyboard
        if ( n==2 && ~isreal(alpha2(index))) || (alpha2(index)<0 )
%             keyboard
            if abs(RAAN_final - RAAN_ToF + f_direction *  2 * (n-1) * pi) > abs(DeltaRAAN_max) || ...
                (alpha2(index)<0 )
                warning('It is not possible to realise this transfer in the given time')
                DeltaV2nd(index) = NaN;
                exitflag = 0;
%                                 keyboard
                continue
            end
        end
        
    end
    
% keyboard
    
    if exitflag == 1

        % DeltaV for the transfer
%         DeltaV2nd(index) = sin(i_final) * sqrt(mu/a_final) *...
%             alpha2(index) / sin(alpha2(index)) * abs(RAAN_final - RAAN_ToF);
        
%         DeltaV2nd(index)= 2 * f * alpha2(index) / pi * abs(RAAN_final - RAAN_ToF) * ...
%             (2 * f * sin(alpha2(index)) / (pi * sin(i_final) ) * sqrt(a_final / mu) -...
%             RAAN_dot)^(-1)
        DeltaV2nd(index)= 2 * f * alpha2(index) / pi * T2nd(index);
         
%          keyboard
        
        
        if DeltaV2nd(index)<0
            keyboard
        end

    else
        alpha2(index) = NaN;
        
    end
    
    
    
    % Total DeltaV, sum of the DeltaV of the 1st and 2nd phases
    DeltaV_TOT(index) = DeltaV1st(index) + DeltaV2nd(index);
    
end

all_solutions.ToF_1st = T1st;
DeltaV_TOT = DeltaV_TOT;

% Find minimum DeltaV
[DeltaV, ind] = (min(DeltaV_TOT));

% Semiamplitude thrust arc first phase of the minimum DeltaV solution
alpha1_min = alpha1(ind);

% Semiamplitude thrust arc second phase of the minimum DeltaV solution
alpha2_min = alpha2(ind);

% Elevation angle of the first phase
beta_min = betaa(ind);

% Time of flight of the first and second phase
T1st_min = T1st(ind);
T2nd_min = T2nd(ind);

% Collect all solutions
all_solutions.DeltaV_TOT = DeltaV_TOT;
all_solutions.alpha1 = alpha1;
all_solutions.alpha2 = alpha2;
all_solutions.DeltaV1 = DeltaV1st;
all_solutions.DeltaV2 = DeltaV2nd;


% Collect results for best solution
best_solution.DeltaV = DeltaV;
best_solution.alpha1_min = alpha1_min;
best_solution.alpha2_min = alpha2_min;
best_solution.beta = beta_min;
best_solution.T1st_min = T1st_min;
best_solution.T2nd_min = T2nd_min;


