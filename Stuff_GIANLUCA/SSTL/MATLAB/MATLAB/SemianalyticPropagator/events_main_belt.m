function [value,isterminal,direction] = events_main_belt(t, x, condition_stop, condition)



a = x(1);

e = sqrt(x(2)^2 + x(3)^2);

ra = a * (1+e);

if strcmp(condition, 'semimajoraxis')
    error = a - condition_stop;
elseif strcmp(condition, 'eccentricity')
    error = e - condition_stop;
elseif strcmp(condition, 'apogee')
    error = ra - condition_stop;
end


value = error;

isterminal = 1;

direction = 0;



return





