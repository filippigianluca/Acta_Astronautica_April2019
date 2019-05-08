function [M,P_imaging,P_idle,info] = space_payload_old_2(x,ep, par)
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
tt        = x(1)*86400;   % actual time, day -> sec         
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
beta_tt = 360*(tt/T_orbit - floor(tt/T_orbit));

% number of completed loops
N_loops = floor(tt/T_orbit);

% time of the last not completed loop
tt_left = (tt-N_loops*T_orbit);






%% EVALUATE DISTANCE CS - GS
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

% cubesat keplerian parameters depending on time
params_CS_tt = [a e I RAAN peri th];

% ground station keplerian parameters (not depending on time)
% params_GS = [R/2 e I0 RAAN peri 0];

% sun keplerian parameters (not depending on time)
% params_SUN = [150000000 e 0 0 0 0];

% %--------------------------------------------------------------------------
% % CARTESIAN orbit parameters
% %--------------------------------------------------------------------------
% % keplerian --> cartesian for the cubesat (CS)
% cart_par_CS = kep2cart(params_CS_tt, mu);
% X_CS = cart_par_CS(1:3);
% 
% % keplerian --> cartesian for the ground station (GS)
% cart_par_GS = kep2cart(params_GS, mu);
% X_GS = cart_par_GS(1:3);
% 
% % distance betweeen CS and GS depending on time
% d_tt = (sum(X_CS - X_GS)^2)^0.5;
% 
% 
% % maximum distance for communication between GS and CS
% d_max_eps0 = R*((((R+h)/R)^2 - cosd(eps0)^2)^0.5 - sind(eps0));
% 
% % maximum distance for communication between GS and CS + take into account
% % uncertainty on the inclination
% d_max_eps0_dI = d_max_eps0*cosd(delta_I*90);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eq = kep2eq(params_CS_tt);
% 
% [L_i, L_o, kep, dt, kep2, alpha_Sun, beta_Sun, phi_Sun] = sol_Eclipse_Eq(Eq, tt, mu, R);
% 
% L=params_CS_tt(6);
% 
% 			dLi = mod(L-L_i, 2*pi);
% 			dLe = mod(L_o-L_i, 2*pi);
% 
% info.L =L;




%--------------------------------------------------------------------------
% propov
%--------------------------------------------------------------------------

% Orbit keplerian elements:
kepel = [a*1000 e I RAAN peri th0];   % see function delkep

% % Orbit state vector:
% stat = kepel_statvec(kepel);
    
% Compute the variations in keplerian elements due to the Earth oblateness
delk = delkep(kepel);

% Orbit propagation
kep2 = kepel + delk*tt;

% To convert from keplerian elements to state vector (if needed)
% stat = kepel_statvec(kep2);






% % Ephemerides date in Modified Julian date:
% year = 2009;
% mjd = djm(13, 4, year);     % Format (day, month, year)
% mjdo = djm(1, 1, year);     % modified julian date of 1/1/year
% mjd1 = djm(1, 1, year+1);   % modified julian date of 1/1/(year+1)
% year_frac = year + (mjd - mjdo)/(mjd1 - mjdo);  % year and fraction
% 
% % Ephemerides time:
% dfra = time_to_dayf (10, 20, 0);    % UTC time in (hour, minute, sec)
% 
%     % Earth's magnetic field
%     tsgr = gst(mjd, dfra + x(1));
%     geoc = inertial_to_terrestrial(tsgr, stat);
%     sphe = rectangular_to_spherical(geoc);
%     alt = sphe(3)/1000;
%     elong = sphe(1);
%     colat = pi/2 - sphe(2);
%     earth_field = 1.e-9*igrf_field (year_frac, alt, colat, elong);  % no sistema NED
%     earth_mag   = rotmay(sphe(2) + pi/2)*rotmaz(-elong)*earth_field'; % para o sistema terrestre
%     earth_mag   = [earth_mag', 0, 0, 0];
% earth_iner  = terrestrial_to_inertial(tsgr, earth_mag); % para o sistema inercial
% sun_sat = sun(mjd, dfra + x(1));
% 
% 
% info.shadow = earth_shadow (kep2, sun_sat);



% L = sum(kep2(4:6));
% info.L = mod(L, 2*pi);


%--------------------------------------------------------------------------
% carlos
% kep2(4:6) = mod(kep2(4:6)./2./pi.*360, 2*pi);
% kep2(1) = kep2(1)/1000;




% params_CS_tt(end) = th0;
Eq = kep2eq(kep2);

% evaluate entring and exit in the shadow
[L_i, L_o, kep, dt, kep2_, alpha_Sun, beta_Sun, phi_Sun] = sol_Eclipse_Eq(Eq, tt_MJD, mu, R);


Li = mod(L_i, 2*pi);
Lo = mod(L_o, 2*pi);

% L = Eq(6);
% dLi = mod(L-L_i, 2*pi);
% dLe = mod(L_o-L_i, 2*pi);

% if 0<= dLi && dLi <= dLe
%     info.L = dLi;
% else
%     info.L = 0;
% end
% info.dt = dt;

% info.shadow = L_o-L_i;
%--------------------------------------------------------------------------

Eclipse = (Lo - Li)/(2*pi); % (% of the day)





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
Tac = 2*beta0_dI/360*T_orbit*N_loops + T_ac_more;  %(sec)


   



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
info.frame_rate = max(info.frame_rate, v_CS/height_coverege);



% vector of number pictures for each loop
N_picture = [];
if floor(tt/T_orbit)>=1
    for i = 1: floor(tt/T_orbit)
        N_picture =  [N_picture floor(info.frame_rate*T_orbit*(1-Eclipse))]; % only with light
    end
end
N_picture = [N_picture floor(info.frame_rate*T_ac_more)];


% N_foto2obdh = info.frame_rate*N_loop + info.frame_rate*more_time*(more_time>T_orbit-Tac);
N_foto_tot  = sum(N_picture); %

% number of images to be stored (conservative quantity)
N_foto2obdh2store = max(N_picture);

% data volume from the camera to store in the OBDH
Data_Volume_to_store = info.image_size*bit_depth/8/2^30*N_foto2obdh2store;  % (Giga bytes)

% total data volume in the mission from the camera 
Data_Volume_tot = info.image_size*bit_depth/8/2^30*N_foto_tot;  % (Giga bytes)




%--------------------------------------------------------------------------
info.N_foto2obdh = N_foto2obdh2store;
info.Data_Volume_obdh_store = info.image_size*bit_depth/8/2^30*max(N_picture); % gigabytes
info.N_foto_tot  = N_foto_tot;

info.Tac_tot         = Tac;
info.Data_Volume            = info.image_size*bit_depth*N_foto_tot;  % (bites)
info.Data_Volume_to_antenna = info.image_size*bit_depth*N_foto_tot;  % (bites)

info.V_tot      = Data_Volume_tot;                % (Giga Bytes)
info.V_to_store = Data_Volume_to_store;           % (Giga Bytes)

info.coverage = coverage;
info.beta0_dI = beta0_dI;

end