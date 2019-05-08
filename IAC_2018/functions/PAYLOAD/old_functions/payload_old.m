function [M, P, info] = payload_old(d, u, par)

Altitude = 10;       % fixed
sensor_area = 10;    % uncertain
% Focal_length = 10;   % uncertain
Circonference_Hearth = 40000; % [km] constant
camera_resolution    = 10;% design
time_mission = 10;   % fixed


% failure parameters
if camera_resolution <= 0.33
    % 1024X1024
    Pixel_size = 4.88;   %[micrometri]
    Focal_length = 19.5; % [mm]
elseif camera_resolution <= 0.66
    % 2048X2048
    Pixel_size = 2.44;   %[micrometri]
    Focal_length = 9.8;  % [mm]  
elseif camera_resolution <= 1
    % 2048X2048
    Pixel_size = 1.22;   %[micrometri]
    Focal_length = 4.9;  % [mm]  
end


Size_picture = Altitude*sensor_area/Focal_length;
length_pictured = Size_picture^0.5;
N_pixel = sensor_area/Pixel_size;


frequence_pictures = length_pictured/Circonference_Hearth;


N_pictures = frequence_pictures*time_mission;

P = 0.5 + picture2power(N_pictures) + time_mission + N_pixel; 
M = 0.5 + picture2mass(N_pictures)  + time_mission + N_pixel; 
end

function F = picture2mass(N_pictures)

F = N_pictures/100;
end

function F = picture2power(N_pictures)

F = N_pictures/100;
end