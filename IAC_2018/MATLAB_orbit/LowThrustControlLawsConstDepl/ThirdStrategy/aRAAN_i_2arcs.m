function [DeltaV, alpha1, ToF1, alpha2, ToF2] = aRAAN_i_2arcs(initial_orbit, final_orbit, ToF_total, ...
                                    f, constants)

% Variation of semimajor axis, inclination and right ascension in given 
% total time of fligth
% Input: initial_orbit -> structure with the initial orbital elements
%        final_orbit   -> structure with the final orbital elements
%        ToF_total     -> total time of fligth
%        f             -> low-thrust acceleration
%        constants     -> structure with constants (J2, gravitational par.)
%        plot_flag     ->

% Output: DeltaV -> cost of the transfer
%         alpha1 -> semiamplitude of the thrust arcs for the first phase
%         ToF1   -> time of flight of the first phase
%         alpha2 -> semiamplitude of the thrust arcs for the second phase
%         ToF2   -> time of flight of the second phase

% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

%% Initialisation

a_initial    = initial_orbit.a;
i_initial    = initial_orbit.incl;
RAAN_initial = initial_orbit.RAAN;

a_final    = final_orbit.a;
i_final    = final_orbit.incl;
RAAN_final = final_orbit.RAAN;

R_Earth = constants.R_Earth;
J2      = constants.J2;
mu      = constants.mu;
                                
                                
%% Find alpha1 and alpha2

n = 0;
for k = 1 : 10
    % --------------------------------------------------------------------
    % Semiamplitude thrust arc second phase
    % --------------------------------------------------------------------
    num2 = 3/4 * pi * mu * J2 * R_Earth^2 * (sin(i_initial) - sin(i_final)) / (f * a_final^4) + ...
        -3/32 * pi * mu * J2 * R_Earth^2 * cos(i_initial) *  (1/a_final^4 - 1/a_initial^4) * ...
        1 / (sqrt(mu/a_initial) - sqrt(mu/a_final)) * (i_final - i_initial) / f * sqrt(mu/a_final);
    
    
    den2 = RAAN_final - RAAN_initial - 2 * n * pi - 3/16 * mu * J2 * R_Earth^2 * cos(i_initial) * ...
        (1/a_final^4 - 1/a_initial^4) * ToF_total * 1 / (sqrt(mu/a_initial) - sqrt(mu/a_final));
    
    alpha2 = asin(num2/den2);
    
    
    % --------------------------------------------------------------------
    % Semiamplitude thrust arc first phase
    % --------------------------------------------------------------------
    num1 = pi/(2*f) * (sqrt(mu/a_initial) - sqrt(mu/a_final));
    den1 = ToF_total - (i_final - i_initial) / (2 * f / pi * sqrt(a_final / mu) * sin(alpha2));
    
    alpha1 = num1 / den1;
    
    % Time of flight second and first phases
    ToF2 = (i_final - i_initial) / (2 * f / pi * sqrt(a_final / mu) * sin(alpha2));
    
    ToF1 = pi/(2 * f * alpha1) * (sqrt(mu/a_initial) - sqrt(mu/a_final));
    
    n = n+1;
    
    if ToF1>0 && ToF2>0 && alpha1<=pi/2 && isreal(alpha2)
        exitflag = 1;
        break
        
    end
        
end

if k == 10 && ~exist('exitflag','var')
    exitflag = 0;
end

if exitflag == 0
    warning('Transfer not possible in given time of flight?')
    DeltaV = NaN;
    return
end


% Time of flight first phase
ToF1 = pi/(2 * f * alpha1) * (sqrt(mu/a_initial) - sqrt(mu/a_final));

% Cost of the trasnfer first and second phase
DeltaV1 = sqrt(mu/a_initial) - sqrt(mu/a_final);
DeltaV2 = alpha2 / sin(alpha2) * sqrt(mu/a_final) * (i_final - i_initial);

% Total DEltaV
DeltaV = DeltaV1 + DeltaV2;