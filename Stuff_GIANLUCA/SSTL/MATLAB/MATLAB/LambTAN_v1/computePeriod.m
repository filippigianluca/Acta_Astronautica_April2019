function [period] = computePeriod(keporb, mu)
% Compute the Period of a keplerian orbit

    % Convert the Semimajor from AU to KM
    a = keporb.a * AstroConstants.Astronomical_Unit;

    % Compute the mean motion [1/sec]
    n = sqrt(mu / a^3);

    % Compute the period [sec]
    period = 2*pi/n;

end
