function R = rotz(z)

% Rotation matrix around the z-axis
% Input angle in degrees

R = [cosd(z)    -sind(z)    0;
     sind(z)     cosd(z)    0;
     0           0          1];

end