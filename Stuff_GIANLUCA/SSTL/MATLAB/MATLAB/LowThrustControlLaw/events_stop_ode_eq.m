function [value,isterminal,direction] = events_stop_ode_eq(t, x, element_target, element_value)



switch element_target
    case 'a'
        error = element_value - x(1); 
    case 'eccentricity'
        error = element_value - sqrt(x(2)^2 + x(3)^2); 
    case 'omega'
        RAAN = atan2(x(4),x(5));
        error = tan(element_value + RAAN) - x(2)/x(3);
    case 'inclination'
        error = tan(element_value/2) - sqrt(x(4)^2 + x(5)^2);
    case 'Omega'
        error = tan(element_value) - x(4)/x(5);
end


value = error;

isterminal = 1;

direction = 0;



return





