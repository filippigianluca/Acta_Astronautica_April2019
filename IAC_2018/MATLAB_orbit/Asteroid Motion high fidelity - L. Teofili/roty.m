function R = roty(y)

% Rotation matrix around the y-axis
% Input angle in degrees

R = [cosd(y)    0     sind(y);
     0          1     0;
    -sind(y)    0     cosd(y)];

end