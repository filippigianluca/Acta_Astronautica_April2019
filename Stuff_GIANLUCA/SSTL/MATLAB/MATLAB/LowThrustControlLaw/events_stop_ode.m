function [value,isterminal,direction] = events_stop_ode(t, x, element_target, element_value)

% Marilena Di Carlo, 2015
% marilena.di-carlo@strath.ac.uk


switch element_target
    case 'a'
        error = element_value - x(1); 
    case 'eccentricity'
        error = element_value - x(2); 
    case 'omega'
        error = element_value - x(5);
    case 'inclination'
        error = element_value - x(3);
    case 'Omega'
        error = element_value - x(4);
end


value = error;

isterminal = 1;

direction = 0;



return





