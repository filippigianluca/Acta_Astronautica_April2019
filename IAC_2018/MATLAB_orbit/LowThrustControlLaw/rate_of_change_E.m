function dx_dt = rate_of_change_E(t, x,  f, beta, constants, mode)

% Equations for the rate of change of the orbital elements with time.
% Input: t-> time
%        x -> vector of orbital elements (a, e, i, RAAN, omega, E)
%        f -> low thrust acceleration
%        beta -> elevation angle
%        constants ->
%        mode -> string for the definition of the time of control
% Output: dx_dt ->
% 
% Referneces:
% Pollard, Simplified Analysis of low-thrust orbital manuevers
% Burt, On space moneuvres with continuous thrust

% Marilena Di Carlo, 2015
% marilena.di-carlo@strath.ac.uk

%%%%%%%%
% f = engine.f;
% m = x(7);
% f = engine.T / m*1e-3;
%%%%%%%%%%%%%%%%%%


% Semimajor axis
a = x(1);

% Eccentricity
e = x(2);

% Inclination
i = x(3);

% RAAN
Omega = x(4);

% omega
omega = x(5);

% Eccentric anomaly
E = mod(x(6), 2*pi);

% True anomaly
cos_theta = (cos(E) - e) / (1 - e * cos(E));
sin_theta = sin(E) * sqrt(1-e^2) / (1 - e * cos(E));
theta = atan2(sin_theta, cos_theta);

% angular momentum
p = a * (1-e^2);
h = sqrt(constants.mu_dim * p);

% Radius
r = p / (1 + e * cos(theta));

dE_dt = constants.mu_dim * sin(theta) / (a * h * sin(E));

% In plane steering
% f1 = -f * (cos(E) - e) / (1 - e * cos(E));
% f2 = f * sqrt(1-e^2) * sin(E) / (1 - e * cos(E));
% f3 = 0;

% % Tangent to orbit
% f1 = f * e * sin(E) / sqrt(1-e^2 * cos(E)^2);
% f2 = f * sqrt(1-e^2)/ sqrt(1-e^2 * cos(E)^2);
% f3 = 0;


switch mode
    
    % Variation of eccentricity with no variation of semimajor axis - Burt
    case 'eccentricity_Delta_a=0'
        f1 = 0;
        f2 = sign(cos(E)) * f ;
        f3 = 0;
        
    % Variation of semimajor axis with no variation of eccentricity - Burt 
    case 'semimajoraxis_Delta_e=0'
        norm_f2 = f / sqrt(1 + 9 * pi * pi * e^2 / (16 * (1-e^2)));
        norm_f1 = 3 * pi * e * norm_f2 / (4 * sqrt(1-e^2));
        f1 = sign(sin(E)) * norm_f1;
        f2 = norm_f2;
        f3 = 0;
    
    %     
    case 'optimum_a'
        if e~=0
        alpha = atan(e * sin(theta) / (1 + e * cos(theta)) );
        else
            alpha = 0;
        end
        f1 = f * sin(alpha);
        f2 = f * cos(alpha);
        f3 = 0;
%         f1 = 0;
%         f2 = f;
%         f3 = 0;
        
    % Variation of eccentricity and inclination - Pollard    
    case 'eccentricity_inclination'
        f1 = f * cos(beta) * sqrt(1-e^2) * sin(E) / (1 - e * cos(E));
        f2 = f * cos(beta) * (cos(E)-e) / (1 - e * cos(E));
        f3 = f * sign(cos(theta)) * sin(beta);
        
    % Maximum variation of eccentricity    
    case 'optimum_eccentricity'
        f1 = f * sin(theta);
        f2 = f * (cos(E) + cos(theta));
        f3 = 0;
        
    % Maximum variation of omega    
    case 'optimum_omega'
        sin_alpha = - cos(theta);
        cos_alpha = (p+r)/p* sin(theta);
%         alpha = atan2(sin_alpha,cos_alpha);
        alpha = atan2(-cos(theta) * (1+e * cos(theta)), (sin(theta) * (2+ e * cos(theta))));
        f1 = f * sin(alpha);
        f2 = f * cos(alpha);
        f3 = 0;
    
    % Maximum variation of omega - Ruggiero    
    case 'optimum_omega_Ruggiero'
        sin_alpha = (1 + e * cos(theta)) * cos(theta);
        cos_alpha = (2 + e * cos(theta)) * sin(theta);
        alpha = atan2(sin_alpha, cos_alpha);
        
        sin_beta = e * sin(omega+theta) * cos(t);
        cos_beta = sin(i) * (sin(alpha-theta) * (1+e*cos(theta)) - cos(alpha)*sin(theta));
        beta = atan2(sin_beta, cos_beta);
        
        f1 = -f * sin(alpha) * cos(beta);
        f2 = -f * cos(alpha) * cos(beta);
        f3 = -f * sin(beta);

    case 'omega'
        f1 = -f * (cos(E) - e) / (1 - e * cos(E));
        f2 = f * sqrt(1-e^2) * sin(E) / (1 - e * cos(E));
        f3 = 0;
        
    case 'inclination'
        f1 = 0;
        f2 = 0;
        f3 =  sign(cos(omega+theta)) * f;
        
    case 'Omega'
        f1 = 0;
        f2 = 0;
        f3 =  sign(sin(omega+theta)) * f;
        
    % Burt 
    case 'omega_Delta_a_Delta_e=0'
        f1 = f;
        f2 = 0;
        f3 = 0;
        
    % Burt 
    case 'omega_Delta_a_Delta_e=0_2'
        f1 = 0;
        f2 = sign(sin(E)) * f;
        f3 = 0;
end


% Rates of change
da_dE = 2 * a^3 / constants.mu_dim * (f1 * e * sin(E) + f2 * sqrt(1-e^2));
da_dt = da_dE * dE_dt;

de_dE = a^2 / constants.mu_dim * (f1 * (1-e^2) * sin(E) + f2 * sqrt(1-e^2) * (2 * cos(E) - e - e * cos(E)^2));
de_dt = de_dE * dE_dt;

di_dE = a^2 / constants.mu_dim * f3 * (1 - e * cos(E)) * (cos(omega) * (cos(E) - e) / sqrt(1-e^2) + ...
                                                       - sin(E) * sin(omega));
di_dt = di_dE * dE_dt;        

dOmega_dE = a^2 / constants.mu_dim * f3 * (1 - e * cos(E)) / sin(i) * (sin(omega) * (cos(E) - e) / sqrt(1-e^2) + ...
                                                       + sin(E) * cos(omega));
dOmega_dt = dOmega_dE * dE_dt;

domega_dE = a^2 / constants.mu_dim * (1/e * (-f1 * sqrt(1-e^2) * (cos(E) - e) + f2 * (2-e^2 - e *cos(E))*sin(E)) + ...
                                         -f3 * (1 - e * cos(E)) * 1 / tan(i) * ...
                                         (sin(omega) * (cos(E) - e) / sqrt(1-e^2) + sin(E) * cos(omega)));
            
domega_dt = domega_dE * dE_dt;


%%%%%%%%        
% dm_dE = - abs(engine.T)/((constants.g0 ) * engine.Isp) * a * h/constants.mu * (1- e * cos(E))/sqrt(1-e^2);
% dt_dE = a * h/constants.mu * (1- e * cos(E))/sqrt(1-e^2); 
%%%%%%%

dx_dt = [da_dt; de_dt; di_dt; dOmega_dt; domega_dt; dE_dt];

end