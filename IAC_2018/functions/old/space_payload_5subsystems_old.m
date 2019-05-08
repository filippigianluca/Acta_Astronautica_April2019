function [M,P_imaging,P_idle,info] = space_payload_5subsystems_old(x,ep, par)
%
% (1) http://propagation.ece.gatech.edu/ECE6390/project/Sum2015/team3/Imaging.html
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





%% CONSTANTS
mu = astro_constants(13); % G*M
R  = astro_constants(23); % earth radius

I0 = par.I0;              % nominal inclination, 0
MJD = par.MJD;            % initial time in MJD, 6.939792e+03 -1


%% DESIGN
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


%% EPISTEMIC
h       = ep(1);      % altitude
eps0    = ep(2);      % elevation angle 
delta_I = ep(3);      % delta inclination orbit









%% CIRCULAR ORBITS (LEO)
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






%--------------------------------------------------------------------------
% KEPLERIAN orbit parameters
%--------------------------------------------------------------------------
a =   (R+h);
e =    0.02;
I =    I0 + delta_I;
RAAN = 0;
peri = 0;
th =   beta_tt;
th0 = 0;




%--------------------------------------------------------------------------
% propov
%--------------------------------------------------------------------------

% Orbit keplerian elements:
kepel = [a*1000 e I RAAN peri th0];   % see function delkep

% % Orbit state vector:
% stat = kepel_statvec(kepel);
    
% Compute the variations in keplerian elements due to the Earth oblateness
delk = delkep(kepel);





%--------------------------------------------------------------------------
% sample the total period of time to evaluate the eclipse time and then the
% volume data got in one orbit



for  k = 1:tt_day/150:tt_day
    
    % Orbit propagation
    kep2 = kepel + delk*tt_day*86400;

    % To convert from keplerian elements to state vector (if needed)
    % stat = kepel_statvec(kep2);

    % params_CS_tt(end) = th0;
    kep2(1) = kep2(1)/1000;
    Eq = kep2eq(kep2);

    % evaluate entring and exit in the shadow
    [L_i, L_o, kep, dt, kep2_, alpha_Sun, beta_Sun, phi_Sun] = sol_Eclipse_Eq(Eq, tt_MJD, mu, R);


    Li = mod(L_i, 2*pi);
    Lo = mod(L_o, 2*pi);

    Eclipse = (Lo - Li)/(2*pi);  % (\% of the day)
    if isnan(Li) || isnan(Lo)
        Eclipse = 0;
    end

    day_time(k) = T_orbit - Eclipse;
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FRACTION OF ACCESS TIME FOR COMMUNICATION FOR EACH LOOP
% nadir angle
alpha0 = asind(R/(R+h)*cosd(eps0));

% central angle
beta0 = 90 - eps0 - alpha0;

% central angle + take into account uncertainty on the inclination
beta0_dI = beta0*cosd(delta_I*90);

% coverage as percentage of earth's surface
coverage = 2*pi*R^2*(1-cosd(beta0))/(4*pi*R^2);


% This term takes into account the not-completed loop
if tt_left/T_orbit >= beta0_dI/360 && tt_left/T_orbit < (360-beta0_dI)/360
    Tac_rate = beta0_dI/360;
elseif tt_left/T_orbit >= (360-beta0_dI)/360
    Tac_rate = beta0_dI/360 + (360*tt_left/T_orbit-(360-beta0_dI))/360;
elseif tt_left/T_orbit < beta0_dI/360
    Tac_rate = (360*tt_left/T_orbit)/360;
end

T_ac_more = Tac_rate*T_orbit; 


% access time: time when the CubeSat is visible from the ground station.
% The first term is the fraction of the central angle and 360 degrees
% multiplied by the number of loops completed. The second term takes into
% account the not-completed loop
Tac = abs(2*beta0_dI/360*T_orbit*N_loops) + T_ac_more;  %(sec)


   



%--------------------------------------------------------------------------
% output

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



% vector of number pictures for each loop
N_images = [];
if tt_sec > T_orbit
%     for i = 1: floor(tt/T_orbit)
        N_images =  [N_images ones(1,floor(tt_sec/T_orbit))*floor(info.frame_rate*T_orbit*(1-Eclipse))]; % images only with light
%     end
end
N_images = [N_images floor(info.frame_rate*T_ac_more)];


% N_foto2obdh = info.frame_rate*N_loop + info.frame_rate*more_time*(more_time>T_orbit-Tac);
N_images_tot  = sum(N_images); %

% number of images to be stored (conservative quantity: I design the compression with the worst case, the biggest number of images in the loops)
N_images_to_compress_and_store = max(N_images);

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

info.coverage = coverage;
info.beta0_dI = beta0_dI;

end