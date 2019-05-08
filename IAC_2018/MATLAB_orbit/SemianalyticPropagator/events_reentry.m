% Function to stop integration when a certain condition on the eccentricity
% or the perigee radius is reached
% Input: t-> time
%        x -> equinoctial state vector
%        condition_reentry -> specify value to stop the integration
%        condition         -> string defining which parameters define when
%        the integration should stop

% Marilena Di Carlo
% marilena.di-carlo@strath.ac.uk


function [value, isterminal, direction] = events_reentry(t, x, condition_reentry,condition)


% Semimajor axis
a = x(1);

% Eccentricity
e = sqrt(x(2)^2 + x(3)^2);

% Perigee radius
rp = a * (1-e);

% Re-entry??
if strcmp(condition, 'altitude')
    error = rp - condition_reentry;
elseif strcmp(condition, 'eccentricity')
    error = e - condition_reentry;
end


value = error;

isterminal = 1;

direction = 0;

return





