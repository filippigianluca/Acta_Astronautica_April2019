function R = rotx(x)

% Rotation matrix around the x-axis
% Input angle in degrees

R = [1      0           0;
     0      cosd(x)    -sind(x);
     0      sind(x)     cosd(x)];
 
end