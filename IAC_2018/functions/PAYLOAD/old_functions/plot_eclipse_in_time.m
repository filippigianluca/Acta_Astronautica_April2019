
clear all; close all; clc


flag = 1;

%% CONSTANT
mu = astro_constants(13); % G*M
R  = astro_constants(23); % earth radius
I0 = 0;                  % nominal inclination

MJD = 6.939792e+03 -1;


%% DESIGN
tt        = 0;        % time (gg)
tt_JD = MJD + tt/86400;
camera = 'meisei_SXGA';
bit_depth = 3;

%% EPISTEMIC
h    = 800;      % altitude, epistemic (km)
eps0 = 0;      % elevation angle
delta_I = 0;   % delta inclination orbit (%)





%--------------------------------------------------------------------------


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
params_GS = [R/2 e I0 RAAN peri 0];

% sun keplerian parameters (not depending on time)
params_SUN = [150000000 e 0 0 0 0];
%--------------------------------------------------------------------------
% CARTESIAN orbit parameters
%--------------------------------------------------------------------------
% keplerian --> cartesian for the cubesat (CS)
cart_par_CS = kep2cart(params_CS_tt, mu);
X_CS = cart_par_CS(1:3);

% keplerian --> cartesian for the ground station (GS)
cart_par_GS = kep2cart(params_GS, mu);
X_GS = cart_par_GS(1:3);

% distance betweeen CS and GS depending on time
d_tt = (sum(X_CS - X_GS)^2)^0.5;


% maximum distance for communication between GS and CS
d_max_eps0 = R*((((R+h)/R)^2 - cosd(eps0)^2)^0.5 - sind(eps0));

% maximum distance for communication between GS and CS + take into account
% uncertainty on the inclination
d_max_eps0_dI = d_max_eps0*cosd(delta_I*90);







if flag == 1
    %--------------------------------------------------------------------------
    % propov
    %--------------------------------------------------------------------------
    
    % Orbit keplerian elements:
    kepel = [a*1000 e I RAAN peri 0];   % see function delkep
    
    % % Orbit state vector:
    % stat = kepel_statvec(kepel);
    
    % Compute the variations in keplerian elements due to the Earth oblateness
    delk = delkep(kepel);
    
    
    % global Number
    
    
    %------------
    step_propagation = 1:T_orbit;
    for aa = step_propagation
        % Number = aa;
        % Orbit propagation
        
%         %---------circular orbit
%         kepel(end) = kepel(end) + aa/T_orbit*2*pi;
%         kepel(1) = kepel(1)/1000;
%         Eq = kep2eq(kepel);
%         tt_JD = tt_JD + aa/(60*60*24);
%         %-----------------------
        
        %---------propagation
        kep2 = kepel + delk*aa;
        kep2(1) = kep2(1)/1000;
        Eq = kep2eq(kep2);
        tt_JD = tt_JD + aa/(60*60*24);
        %------------------------
        [L_i, L_o, kep, dt, kep2_, alpha_Sun, beta_Sun, phi_Sun] = sol_Eclipse_Eq(Eq, tt_JD, mu, R);
        tempo(aa) = dt;
    end

   
%     plot(plot_period, ones(1,plot_period(end)).*T_orbit, 'linewidth', 2)
    hold on
    plot(step_propagation, T_orbit - tempo(step_propagation), 'linewidth', 2)
    %------------
    
    
    
    
    
elseif flag == 2 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    step_propagation = 1:10000;
    for aa = step_propagation
        
        tt        = aa;
        % central angle function of time
        beta_tt = 360*(tt/T_orbit - floor(tt/T_orbit));
        
        % number of completed loops
        N_loops = floor(tt/T_orbit);
        
        % time of the last not completed loop
        tt_left = (tt-N_loops*T_orbit);
        
        
        
        
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
        
        
        
        tempo_plot(aa) = Tac;
    end
    
    plot(tempo_plot)

    
elseif flag == 3
    
    
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


    
    
    
    
    %----------
     step_propagation = 1:10000;
  for aa = step_propagation
        
        tt        = aa;
        % central angle function of time
        beta_tt = 360*(tt/T_orbit - floor(tt/T_orbit));
        
        % number of completed loops
        N_loops = floor(tt/T_orbit);
        
        % time of the last not completed loop
        tt_left = (tt-N_loops*T_orbit);
        
        
        kepel = [a*1000 e I RAAN peri 0];
        delk = delkep(kepel);
        
        
        kep2 = kepel + delk*aa;
        kep2(1) = kep2(1)/1000;
        Eq = kep2eq(kep2);
        tt_JD = tt_JD + T_orbit/(60*60*24);
        [L_i, L_o, kep, dt, kep2_, alpha_Sun, beta_Sun, phi_Sun] = sol_Eclipse_Eq(Eq, tt_JD, mu, R);
        
        Eclipse = dt/T_orbit; % (% of the day)
        
        
        
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
    % total data volume in the mission from the camera 
    Data_Volume_tot = info.image_size*bit_depth/8/2^30*N_foto_tot;  % (Giga bytes)    
    
    
    
    output(aa)         = Data_Volume_tot/Tac;
    
    outputDV(aa)       = Data_Volume_tot;
    outputTac(aa)      = Tac;
    output_eclipse(aa) = Eclipse;
    
    end   
    %----------
    
    
    


    
end


