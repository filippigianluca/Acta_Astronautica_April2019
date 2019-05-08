function [M,P_imaging,P_idle,info] = IAC_payload(x,ep, par)
%
% (1) http://propagation.ece.gatech.edu/ECE6390/project/Sum2015/team3/Imaging.html
%
% the length of the eclipse is here function only of the altitude as in:
% (2) https://www.quora.com/If-satellite-is-revolving-700-km-above-Earth-then-what-will-be-the-eclipse-duration-for-a-satellite-What-happens-to-the-eclipse-duration-if-the-satellite-altitude-is-increased
%
%
% INPUT: - design
%           * x(1): tt, time of the mission    (days);
%           * x(2): depth for pixel definition (integer);
%           * x(3): type of camera;
%        - epistemic
%           * ep(1): h, altitude,           (km);
%           * ep(2): eps0, elevation angle, (deg);
%           * ep(3): percentage of inclination of the orbit (Delta %I)
%
% OUTPUT:   * M:         mass of the camera;
%           * P_imaging: power needed during the light time;
%           * P_idle:    power needed during night;
%           * info:      -




%--------------------------------------------------------------------------
% FIXED

cost_r1 = 0.5;

mu = astro_constants(13); % G*M
R  = astro_constants(23); % earth radius

I0 = par.I0+30;              % nominal inclination, 0
MJD = par.MJD;            % initial time in MJD, 6.939792e+03 -1



lat_GS =   30;            % latitude of the ground station
DELTA_long = 10;          % difference in longitude between orbit pole and ground station
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
% DESIGN

tt_day    = x(1);
tt_sec    = x(1)*86400;   % actual time, day -> sec         
tt_MJD    = MJD + x(1);   % actual time, MJD 
bit_depth = ceil(x(2));    

if x(3)  <= 0.25          
    camera = 'meisei_SXGA';
elseif x(3) <= 0.5
    camera = 'meisei_VGA';
elseif x(3) <= 0.75
    camera = 'ecam_C50';
else 
    camera = 'ecam_dvr4';
end
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
% EPISTEMIC

h            = ep(1);      % altitude
epsilon_min  = ep(2);      % elevation angle 
delta_I      = ep(3);      % delta inclination orbit
%--------------------------------------------------------------------------






% CIRCULAR ORBITS (LEO)
%--------------------------------------------------------------------------
% period of the circolar orbit
T_orbit = 2*pi*((R+h)^3/mu)^0.5;   % (sec)

% veloecity
v_CS = 2*pi*(R+h)*1000/T_orbit;

% central angle function of time
beta_tt = 360*(tt_sec/T_orbit - floor(tt_sec/T_orbit));

% number of completed loops
N_loops = floor(tt_sec/T_orbit);

% time of the last not completed loop
tt_left = (tt_sec-N_loops*T_orbit);

info.N_CL    = tt_day*86400/T_orbit;

% earth angular diameter
EAD = 2*(asind((R/(h+R)))); 

% eclispse in hours
ECLIPSE = EAD/360*T_orbit;


info.t_day   = T_orbit - ECLIPSE; % hours
info.t_night = ECLIPSE;           % hours
%--------------------------------------------------------------------------



% access time
%--------------------------------------------------------------------------
rho = EAD/2;

eta_max = asind(sind(rho)*cosd(epsilon_min));

lambda_max = 90 - epsilon_min - eta_max;

I =    I0 + delta_I;

lat_pole = 90 - I;

lambda_min = asind( sind(lat_pole)*sin(lat_GS) + cosd(lat_pole)*cosd(lat_GS)*cosd(DELTA_long) );


Tac = (T_orbit/180)*acosd(cosd(lambda_max/cosd(lambda_min)));

%--------------------------------------------------------------------------







% choose the payload type
%--------------------------------------------------------------------------
switch camera % from data sheet
    case 'meisei_SXGA'
        M = 1.1;                     % kg
        P_imaging = 4;               % Watt, imaging, day
        P_idle = 0;                  % Watt, idle, night  
        info.image_size = 1280*1024; % pixel 
        info.frame_rate = 6.6;       % sec^-1
    case 'meisei_VGA'
        M = 1.1;                     
        P_imaging = 4;  
        P_idle = 0;             
        info.image_size = 640*480;
        info.frame_rate = 26.6;       
    case 'ecam_C50'
        M = 0.256;                     
        P_imaging = 2.5;  
        P_idle = 1.75;            
        info.image_size = 2592*1944;   
        info.frame_rate = 26.6;       
    case 'ecam_dvr4'
        M = 1.1;                     
        P_imaging = 4; 
        P_idle = 9.75;           
        info.image_size = 1280*1024;  
        info.frame_rate = 26.6;       
end

height_coverege = 23000; % (meters)
info.frame_rate = min(info.frame_rate, v_CS/height_coverege); % from (1)
%--------------------------------------------------------------------------


















%--------------------------------------------------------------------------
% sample the total period of time to evaluate the eclipse time and then the
% volume data got in one orbit


TIME_vector = 0:tt_day/150:tt_day;

dt = 0.1*ones(1,length(TIME_vector));
N_images_vector = floor(info.frame_rate*T_orbit*(1-dt/T_orbit));


R2 = interp1(TIME_vector, N_images_vector, tt_day,'linear','extrap');
R1 = R2*cost_r1;

R1_function = @(x)(R1fun(TIME_vector,N_images_vector,x));
R2_function = @(x)(R2fun(TIME_vector,N_images_vector,x));


global N_images_vector2R_function
N_images_vector2R_function = N_images_vector;

global TIME_vector2R_function
TIME_vector2R_function = TIME_vector;

R2_vector = N_images_vector;
R1_vector = R2_vector.*cost_r1;
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------

% N_foto2obdh = info.frame_rate*N_loop + info.frame_rate*more_time*(more_time>T_orbit-Tac);
N_images_tot  = sum(N_images_vector); %

% number of images to be stored (conservative quantity: I design the compression with the worst case, the biggest number of images in the loops)
N_images_to_compress_and_store = max(N_images_vector);

% data volume from the camera to store in the OBDH
Data_Volume_to_compress_and_store = info.image_size*bit_depth/8/2^30*N_images_to_compress_and_store;  % (Giga bytes)

% total data volume in the mission from the camera 
Data_Volume_tot = info.image_size*bit_depth/8/2^30*N_images_tot;  % (Giga bytes)




%--------------------------------------------------------------------------
info.N_images_to_compress_and_store = N_images_to_compress_and_store;
info.N_images_tot                   = N_images_tot;

info.Data_Volume_tot                   = Data_Volume_tot;                      % (Giga Bytes)
info.Data_Volume_to_compress_and_store = Data_Volume_to_compress_and_store;    % (Giga Bytes)

info.Tac_tot         = Tac;

% info.coverage = coverage;
% info.beta0_dI = beta0_dI;


info.R2 = R2;
info.R1 = R1;

info.R1_vector = R1_vector;
info.R2_vector = R2_vector;

info.R1_function = R1_function;
info.R2_function = R2_function;

end