% function of the cubesat with 3 systems: TTC, AOCS and POWER. I use
% design, uncertain and fixed parameters as in Alicino paper. All the
% parameters not in Alicino paper are fixed with values found in internet.
% POWER system function is modified to be as like in the paper.

function output = IAC2018_only_mass2_for_belief(d, u, par) 




% to change: 
% 1. ep_ttc(4) = 640;        %%%       * ep(4) = range, km  --> should be the altitude
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAYLOAD

% design
x_payload(1) = par.fix.time;   % TIME, fixed (days)
x_payload(2) = d(10);   % bit_depth 
x_payload(3) = d(11);   % camera_catalogs 0.4

% epistemic
ep_payload(1) = u(17);      % h, altitude, epistemic (km)  [600 800][800 1000]
ep_payload(2) = u(18);      % elevation angle (deg)        [0 10][10 20]
ep_payload(3) = u(19);      % delta inclination orbit (%)  [0 10][10 20]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBDH

% design
x_obdh(1) = par.fix.time;  % TIME 
x_obdh(2) = d(12);  % type  [0 1]

% epistemic
% ep_obdh(1); % N_foto2obdh, maximum number of pictures to compress and store in all he loops
% ep_obdh(2); % N_foto_tot,
% ep_obdh(3); % Data_Volume_obdh_store
ep_obdh(5) = u(20); % \% mass   [0 10][10 20]
ep_obdh(6) = u(21); % \% power  [0 10][10 20]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AOCS

% design
x_aocs(11) =  d(1);    %       * x(11) = slew angle, deg                                             (ALICINO design, [10 60])
x_aocs(12) =  d(2);    %       * x(12) = time for slew maneuvre, s                                   (ALICINO design, [30 90])

% uncertain
x_aocs(2) = u(1);      %       * x(2) = offset between centre of gravity and centre of pressure, m   (ALICINO uncertain, [0.005 0.01][0.01 0.02])
x_aocs(3) = u(2);      %       * x(3) = area normal to velocity vector, m^2                          (ALICINO uncertain, [0.034 0.0885][0.0885 0.15])
x_aocs(5) = u(3);      %       * x(5) = reflectance factor                                           (ALICINO uncertain, [0.5 0.6][0.6 0.7])
x_aocs(6) = u(4);      %       * x(6) = spacecraft residual dipole, A*m^2                            (ALICINO uncertain, [0.0005 0.001][0.001 0.0015])
x_aocs(7) = u(5);      %       * x(7) = drag coefficient                                             (ALICINO uncertain, [2 2.2][2.2 2.5])
ep_aocs(14) = u(6);    %       * ep(14) = Delta Inertia                                              (ALICINO uncertain, [-10 5][5 10])  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  TTC 

% design
x_ttc(1) = d(3);       %       * x(1) = frequency, GHz                           (ALICINO design       [7 10])
x_ttc(2) = d(4);       %       * x(2) = modulation,                              (ALICINO design       [0 1])
x_ttc(9) = d(5);       %       * x(9) = amplifier type,                          (ALICINO design       [0 1])

% uncertain
x_ttc(3) = u(7);       %       * x(3) = antenna efficiency                       (ALICINO uncertainty  [0.6 0.8] [0.8 0.9])
x_ttc(4) = u(8);       %       * x(4) = antenna gain, dB                         (ALICINO uncertainty  [1 3][3 5] trasmit AG)
x_ttc(5) = u(9);       %       * x(5) = onboard loss, dB                         (ALICINO uncertainty  [0.1 0.5][0.5 1])
x_ttc(6) = u(10);      %       * x(6) = other unmodelled losses, dB              (ALICINO uncertainty  [0.5 1.5][1.5 2])
x_ttc(7) = u(11);      %       * x(7) = mass of distribution network, kg         (ALICINO uncertainty  [0.1 0.3][0.2 0.5])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  POWER

% design
x_power(1) =  d(6);    %       * x(1) = type of solar cell                                                       (design, [0 1])
x_power(2) =  d(7);    %       * x(2) = required bus voltage, V [0,100]                                  (ALICINO design, [3 5]) 
x_power(3) =  d(8);    %       * x(3) = eps configuration (det, mppt) [0,0.5) = det; [0.5,1] = mppt      (ALICINO design, [0 1])
x_power(8) =  d(9);    %       * x(8) = allowable voltage drop,                                          (ALICINO design, [1 3])              

% uncertain
x_power(4) =  u(12);   %       * x(4)   = cell packing efficiency (or assembly factor) [0,1]     (ALICINO uncertain, [0.8 0.85][0.85 0.9])   
ep_power(10) = u(14);  %       * ep(10) = uncertainty on power requirements, %                   (ALICINO uncertain, [0 10][10 20])

%-- add to the function "space_power_Alicino_paper":
ep_power(13) = u(15);  %        uncertain         Delta rho_sa [% kg/m^2] density solar array, [0 15][15 30]
ep_power(14) = u(16);  %        uncertain         Delta D_cell [%] degradation solar array,    [0 50][50 100]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIXED for AOCS 
%--------------------------------------------------------------------------
x_aocs(1) = 30;         %%%       * x(1) = pointing accuracy, deg                                                (fix, 30)   
x_aocs(4) = x_aocs(3);     %       * x(4) = area normal to sun line, m^2                                 (ALICINO uncertain, [0.034 0.0885][0.0885 0.15])

x_aocs(8) = 0.1;        %%%       * x(8) = actuator type: [0 0.5) = magtorq, [0.5 1] = thruster                    (magnetotorque or thruster, fix, 0.1)          
x_aocs(9) =  200;       %%%       * x(9) = specific impulse, s                                                   (only for thruster, fix 200) 
x_aocs(10) =  50;       %%%       * x(10) = burn time for momentum dumping, s                                    (only for thruster, fix 50) 

x_aocs(13) =  0;        %%%       * x(13) = initial tumbling rate, rad/s                                         (fix, 10)
x_aocs(14) =  0;        %%%       * x(14) = spin angle, deg                                                      (fix, 100)

ep_aocs(1) =  5850;     %%%       * ep(1) = orbit period, s                                             (ALICINO fixed, 5850)           
ep_aocs(2) =  3.1e-5;   %%%       * ep(2) = mean planet's magnetic field stregth, T                     (ALICINO fixed, 3.1e-5)
ep_aocs(3) =  1420;     %%%       * ep(3) = mean incident solar radiation, W/m^2                        (ALICINO fixed, 1420)
ep_aocs(4) =  1/2*1e-5*7.54^2;%%%       * ep(4) = mean dynamic pressure, Pa                                      (1/2 * rho * v^2 = 1/2 * 1e-5  * 7.54^2)
ep_aocs(5) = 1e-5;      %%%       * ep(5) = mean gravitational field strength, deg                               (fixed)
ep_aocs(6) =  5;        %%%       * ep(6) = angle between spacecraft z axis and nadir vector, deg       (ALICINO fixed, 5)
ep_aocs(7) =  0.0417;   %%%       * ep(7) = spacecraft moment of inertia about x axis, kg*m^2           (ALICINO fixed, 0.0417)
ep_aocs(8) =  0.1083;   %%%       * ep(8) = spacecraft moment of inertia about y axis, kg*m^2           (ALICINO fixed, 0.1083)
% ep_aocs(9) =  0.1417;   %%%       * ep(9) = spacecraft moment of inertia
% about z axis, kg*m^2           (ALICINO fixed, 0.1417) IS NOW FUNCTION OF
% THE PAYLOAD
ep_aocs(10) =  1;       %%%       * ep(10) = thruster moment arm, m                                              (only for thruster, fix 1) 
ep_aocs(11) =  1;       %%%       * ep(11) = number of thrusters                                                 (only for thruster, fix 1) 
ep_aocs(12) =  1;       %%%       * ep(12) = number of dumping maneuvres during S/C lifetime                     (fix 1) 
ep_aocs(13) =  1;       %%%       * ep(13) = number of slew maneuvres during S/C lifetime                        (fix 1) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIXED for TTC 
%--------------------------------------------------------------------------
x_ttc(8) = 0.5;         %%%       * x(8) = modulation index, rad

ep_ttc(1) = 1e-5;       %%%       * ep(1) = Bit Error Rate                          (ALICINO fix, 1e-5)
% ep_ttc(2) = 1e6;        %%%       * ep(2) = data volume, bits
% (ALICINO fix, 1e6) FROM PAYLOAD
ep_ttc(3) = 30;         %%%       * ep(3) = ground station G/T, dB                  (ALICINO fix, 30)
ep_ttc(4) = 640;        %%%       * ep(4) = range, km                               (ALICINO fix, 640)
ep_ttc(5) = 10;         %%%       * ep(5) = elevation angle, deg                    (ALICINO fix, 10)
ep_ttc(6) = 5;          %%%       * ep(6) = pointing accuracy, deg                  (ALICINO fix, 5)
c = 299792458;
ep_ttc(7) = (c/(1e9*x_ttc(1)))/pi*(10e6/0.55)^0.5;     %%%       * ep(7) = ground station antenna diameter, m      (lambda/pi * (10e6/0.55)^0.5; lambda = c/(1e9*f))
% ep_ttc(8) = 10;         %%%       * ep(8) = ground station access time, min         (ALICINO fix, 10)
ep_ttc(9) = 0;          %%%       * ep(9) = link margin, dB                         (margin on the P, i can fix to 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIXED for POWER 
%--------------------------------------------------------------------------
x_power(6) =   1;        %%%       * x(6) = number of loads (to determine LCLs in PCDU) [1,50]
x_power(7) =   10;       %%%       * x(7) = harness mass factor (as a percentage of eps mass) (harness=cables and buses)

ep_power(5) = 1;         %%%       * ep(5) = distance from Sun, AU                                 (ALICINO fix, 1)
ep_power(6) = 5;         %%%       * ep(6) = worst case angle of incidence, deg                    (ALICINO fix, 5)
ep_power(7) = 0; %u(13);   %       * ep(7)  = uncertainty on array temperature, C                    (ALICINO uncertain, [0 10][10 15])
ep_power(8) = 1;         %%%       * ep(8) = satellite lifetime, yrs                               (ALICINO fix, 1)
ep_power(9) = 1;         %%%       * ep(9) = bus regulation (0: unregulated, 1: regulated)         (ALICINO fix, 1)

ep_power(11) = 1e14;     %%%       * ep(11) = total radiation fluence                                      (fixed)
ep_power(12) = 1;        %%%      * ep(12) = s/c type (1: cubesat, 0: big one)                     (ALICINO fix, 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIXED for PAYLOAD
%--------------------------------------------------------------------------

par_payloads.MJD = 6.939792e+03 -1;
par_payloads.I0  = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% OBJECTIVE FUNCTION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------






%% -------------------------------------------------------------------------
% Payload: camera
%--------------------------------------------------------------------------
[M_payload, P_payload_imaging, P_payload_idle, info_payload] = space_payload_5subsystems(x_payload, ep_payload, par_payloads);
%--------------------------------------------------------------------------



%% -------------------------------------------------------------------------
% On Board Data Handling
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
ep_obdh(1) = info_payload.N_images_to_compress_and_store;    % maximum number of pictures to compress and store in all he loops
ep_obdh(2) = info_payload.N_images_tot;                      % total number of images
ep_obdh(3) = info_payload.Data_Volume_to_compress_and_store; % data volume to store in OBDH (Giga Bytes)
ep_obdh(4) = info_payload.Data_Volume_tot;                   % data volume total            (Giga Bytes)
%--------------------------------------------------------------------------
[M_obdh, P_obdh, info_obdh] =  space_obdh_5subsystems(x_obdh, ep_obdh);
%--------------------------------------------------------------------------
% 
% 
% 
% 
% 
% 
%% ------------------------------------------------------------------------
% Attitude and Orbit Control
%--------------------------------------------------------------------------
ep_aocs(9) =  0.1417 + M_payload*x_aocs(4);   %       * ep(9) = spacecraft moment of inertia about z axis, kg*m^2           (ALICINO fixed, 0.1417)
ep_aocs(14) = u(17); % altitude h
%--------------------------------------------------------------------------
[M_aocs, P_aocs, info_aocs] = space_aocs_5subsystems(x_aocs, ep_aocs);
% 
% 
% 
% 
% 
% 
%% ------------------------------------------------------------------------
% Telecommunication 
%--------------------------------------------------------------------------
ep_ttc(2) = info_obdh.V_compressed_tot*8*2^30;   %       * ep(2) = data volume, Gbytes -> bits             
ep_ttc(8) = info_payload.Tac_tot;                %       * ep(8) = ground station access time, min         

[M_ttc, P_ttc, info_ttc] = space_ttc_5subsystems(x_ttc,ep_ttc);
%--------------------------------------------------------------------------
% 
% 
% 
% 
%% ------------------------------------------------------------------------
% Power
%--------------------------------------------------------------------------
N_CL = info_payload.N_CL;                            % number of charge/discharge
x_power(5) =   -36.76*log(N_CL/207800)/100;          % x(5) = depth of discharge (it is how much the battery has been discharged) [0,1]        (ALICINO fix, DOD = -36.76 ln(CL/207800), CL 14250)
ep_power(3) = info_payload.t_day/par.fix.time;                    % fixed to 1.615 by alicino     %%%       * ep(3) = orbital daylight time, h                              (ALICINO fix, 1.615)
ep_power(4) = info_payload.t_night/par.fix.time;                  % 0.0103;    %%%       * ep(4) = orbital eclipse time, h                               (ALICINO fix, 0.0103)

% exchange function: power
ep_power(1) = 16+P_aocs + P_ttc + P_payload_imaging + P_obdh;    %       * ep(1) = required daylight power, W                      (ALICINO transfer function)
ep_power(2) = 16+P_aocs + P_ttc + P_payload_idle    + P_obdh;    %       * ep(2) = required eclipse power, W                       (ALICINO transfer function)

[M_power, P_power, info_power] = space_power_5subsystems(x_power,ep_power);
%--------------------------------------------------------------------------
% 
% 
% 
% 
% 
% 
% 
% % sum the output from the subsystems
% %--------------------------------------------------------------------------
% M_tot = M_payload + M_obdh + M_aocs + M_ttc + M_power;
% % P_tot = P_payload_imaging + P_payload_idle + P_obdh + P_aocs + P_ttc + P_power;
% % V = info_obdh.V_compressed_tot;
% %--------------------------------------------------------------------------




% % normalisation of Data Volume and Mass
% %--------------------------------------------------------------------------
% V_reference = 1e6;
% M_tot_reference = 50;
% 
% V_normalised = (V_reference - V)/V_reference;
% M_normalised = (M_tot_reference - M_tot)/M_tot_reference;
% %--------------------------------------------------------------------------



%% RES
% R1_function = info_payload.R1_function;
% R2_function = info_payload.R2_function;
% 
% 
% TM = par.fix.time;
% 
% lam0 = 1;
% mu0 = 2;
% [lam,mu] = cubesat_failure_rate(d,u, lam0, mu0);
% lam = lam/365;
% mu = mu/365;
% 
% RES = resilienceSimpleFunction(TM,R2_function,R1_function,lam, mu);
%--------------------------------------------------------------------------

% minimise mass and maximise RES
% F = - RES;
% F = M_tot;


%--------------------------------------------------------------------------
%% output
%--------------------------------------------------------------------------
output = M_payload + M_obdh + M_aocs + M_ttc + M_power;
% output.c = - RES + par.fix.nu;
% output.ceq = [];
%--------------------------------------------------------------------------

return