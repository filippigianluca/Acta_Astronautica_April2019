function [ hh mm ss ] = seconds2HHMMSS( seconds )

    % Compute the Hours
    hh = floor (seconds/3600);

    % Compute the Minutes
    mm = floor ( (seconds - hh * 3600) / 60);

    % Compute the Seconds
    ss = seconds - mm * 60 - hh * 3600;

end

