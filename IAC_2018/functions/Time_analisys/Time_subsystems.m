% time running functions

run_flag = 6;

% 1 -> payload
% 2 -> obdh
% 3 -> aocs
% 4 -> ttc
% 5 -> power
% 6 -> all 5 functions ("IAC2018_only_mass2_for_belief")
% 7 -> Danda function ("resilienceSimpleFunction")



if run_flag == 1
    %======================================================================
    % PAYLOAD
    %======================================================================
    % design
    x_payload(1) = 365;   % TIME, fixed (days)
    x_payload(2) = rand;   % bit_depth
    x_payload(3) = rand;   % camera_catalogs 0.4
    
    % epistemic
    ep_payload(1) = rand;      % h, altitude, epistemic (km)  [600 800][800 1000]
    ep_payload(2) = rand;      % elevation angle (deg)        [0 10][10 20]
    ep_payload(3) = rand;      % delta inclination orbit (%)  [0 10][10 20]
    
    % par
    par_payloads.MJD = 6.939792e+03 -1;
    par_payloads.I0  = 0;
    
    tic
    for i=1:1000
        [M_payload, P_payload_imaging, P_payload_idle, info_payload] = space_payload_5subsystems(x_payload, ep_payload, par_payloads);
    end
    toc
    
    
elseif run_flag == 2
    %======================================================================
    % OBDH
    %======================================================================
    % design
    x_obdh(1) = 365;  % TIME
    x_obdh(2) = rand;  % type  [0 1]
    
    % epistemic
    ep_obdh(5) = rand; % \% mass   [0 10][10 20]
    ep_obdh(6) = rand; % \% power  [0 10][10 20]
    ep_obdh(1) = rand;    % maximum number of pictures to compress and store in all he loops
    ep_obdh(2) = rand;                      % total number of images
    ep_obdh(3) = rand; % data volume to store in OBDH (Giga Bytes)
    ep_obdh(4) = rand;                   % data volume total            (Giga Bytes)
    %----------------------------------------------------------------------
    tic
    for i=1:1000
        [M_obdh, P_obdh, info_obdh] =  space_obdh_5subsystems(x_obdh, ep_obdh);
    end
    toc
    %----------------------------------------------------------------------
    
    
elseif run_flag == 3
    %======================================================================
    % AOCS
    %======================================================================
    % design
    x_aocs(11) =  rand;    %       * x(11) = slew angle, deg                                             (ALICINO design, [10 60])
    x_aocs(12) =  rand;    %       * x(12) = time for slew maneuvre, s                                   (ALICINO design, [30 90])
    
    % uncertain
    x_aocs(2) = rand;      %       * x(2) = offset between centre of gravity and centre of pressure, m   (ALICINO uncertain, [0.005 0.01][0.01 0.02])
    x_aocs(3) = rand;      %       * x(3) = area normal to velocity vector, m^2                          (ALICINO uncertain, [0.034 0.0885][0.0885 0.15])
    x_aocs(5) = rand;      %       * x(5) = reflectance factor                                           (ALICINO uncertain, [0.5 0.6][0.6 0.7])
    x_aocs(6) = rand;      %       * x(6) = spacecraft residual dipole, A*m^2                            (ALICINO uncertain, [0.0005 0.001][0.001 0.0015])
    x_aocs(7) = rand;      %       * x(7) = drag coefficient                                             (ALICINO uncertain, [2 2.2][2.2 2.5])
    ep_aocs(14) = rand;    %       * ep(14) = Delta Inertia                                              (ALICINO uncertain, [-10 5][5 10])
    ep_aocs(9) =  rand;
    
    x_aocs(1) = 30;         %%%       * x(1) = pointing accuracy, deg                                                (fix, 30)
    x_aocs(4) = x_aocs(3);     %       * x(4) = area normal to sun line, m^2                                 (ALICINO uncertain, [0.034 0.0885][0.0885 0.15])
    
    x_aocs(8) = 0.1;        %%%       * x(8) = actuator type: [0 0.5) = magtorq, [0.5 1] = thruster                    (magnetotorque or thruster, fix, 0.1)
    x_aocs(9) =  200;       %%%       * x(9) = specific impulse, s                                                   (only for thruster, fix 200)
    x_aocs(10) =  50;       %%%       * x(10) = burn time for momentum dumping, s                                    (only for thruster, fix 50)
    
    x_aocs(13) =  0;        %%%       * x(13) = initial tumbling rate, rad/s                                         (fix, 10)
    x_aocs(14) =  0;        %%%       * x(14) = spin angle, deg                                                      (fix, 100)
    %----------------------------------------------------------------------
    tic
    for i=1:1000
        [M_aocs, P_aocs, info_aocs] = space_aocs_5subsystems(x_aocs, ep_aocs);
    end
    toc
    %----------------------------------------------------------------------
    
elseif run_flag == 4
    %======================================================================
    % TTC
    %======================================================================
    x_ttc(1) = rand;       %       * x(1) = frequency, GHz                           (ALICINO design       [7 10])
    x_ttc(2) = rand;       %       * x(2) = modulation,                              (ALICINO design       [0 1])
    x_ttc(9) = rand;       %       * x(9) = amplifier type,                          (ALICINO design       [0 1])
    
    % uncertain
    x_ttc(3) = rand;       %       * x(3) = antenna efficiency                       (ALICINO uncertainty  [0.6 0.8] [0.8 0.9])
    x_ttc(4) = rand;       %       * x(4) = antenna gain, dB                         (ALICINO uncertainty  [1 3][3 5] trasmit AG)
    x_ttc(5) = rand;       %       * x(5) = onboard loss, dB                         (ALICINO uncertainty  [0.1 0.5][0.5 1])
    x_ttc(6) = rand;      %       * x(6) = other unmodelled losses, dB              (ALICINO uncertainty  [0.5 1.5][1.5 2])
    x_ttc(7) = rand;      %       * x(7) = mass of distribution network, kg         (ALICINO uncertainty  [0.1 0.3][0.2 0.5])
    x_ttc(8) = 0.5;
    
    ep_ttc(2) = rand;   %       * ep(2) = data volume, Gbytes -> bits
    ep_ttc(8) = rand;                %       * ep(8) = ground station access time, min
    
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
    %----------------------------------------------------------------------
    tic
    for i=1:1000
        [M_ttc, P_ttc, info_ttc] = space_ttc_5subsystems(x_ttc,ep_ttc);
    end
    toc
    %----------------------------------------------------------------------
    
    
elseif run_flag == 5
    %======================================================================
    % POWER
    %======================================================================
    % design
    x_power(1) =  rand;    %       * x(1) = type of solar cell                                                       (design, [0 1])
    x_power(2) =  rand;    %       * x(2) = required bus voltage, V [0,100]                                  (ALICINO design, [3 5])
    x_power(3) =  rand;    %       * x(3) = eps configuration (det, mppt) [0,0.5) = det; [0.5,1] = mppt      (ALICINO design, [0 1])
    x_power(8) =  rand;    %       * x(8) = allowable voltage drop,                                          (ALICINO design, [1 3])
    
    % uncertain
    x_power(4)   =  rand;   %       * x(4)   = cell packing efficiency (or assembly factor) [0,1]     (ALICINO uncertain, [0.8 0.85][0.85 0.9])
    ep_power(7)  = rand;   %       * ep(7)  = uncertainty on array temperature, C                    (ALICINO uncertain, [0 10][10 15])
    ep_power(10) = rand;  %       * ep(10) = uncertainty on power requirements, %                   (ALICINO uncertain, [0 10][10 20])
    
    %-- add to the function "space_power_Alicino_paper":
    ep_power(13) = rand;  %        uncertain         Delta rho_sa [% kg/m^2] density solar array, [0 15][15 30]
    ep_power(14) = rand;  %        uncertain         Delta D_cell [%] degradation solar array,    [0 50][50 100]
    
    
    x_power(6) =   1;        %%%       * x(6) = number of loads (to determine LCLs in PCDU) [1,50]
    x_power(7) =   10;       %%%       * x(7) = harness mass factor (as a percentage of eps mass) (harness=cables and buses)
    ep_power(5) = 1;         %%%       * ep(5) = distance from Sun, AU                                 (ALICINO fix, 1)
    ep_power(6) = 5;         %%%       * ep(6) = worst case angle of incidence, deg                    (ALICINO fix, 5)
    ep_power(8) = 1;         %%%       * ep(8) = satellite lifetime, yrs                               (ALICINO fix, 1)
    ep_power(9) = 1;         %%%       * ep(9) = bus regulation (0: unregulated, 1: regulated)         (ALICINO fix, 1)
    ep_power(11) = 1e14;     %%%       * ep(11) = total radiation fluence                                      (fixed)
    ep_power(12) = 1;        %%%      * ep(12) = s/c type (1: cubesat, 0: big one)                     (ALICINO fix, 1)
    
    N_CL = rand;                            % number of charge/discharge
    x_power(5) =   -36.76*log(N_CL/207800)/100;          % x(5) = depth of discharge (it is how much the battery has been discharged) [0,1]        (ALICINO fix, DOD = -36.76 ln(CL/207800), CL 14250)
    ep_power(3) = rand;                    % fixed to 1.615 by alicino     %%%       * ep(3) = orbital daylight time, h                              (ALICINO fix, 1.615)
    ep_power(4) = rand;                    % 0.0103;    %%%       * ep(4) = orbital eclipse time, h                               (ALICINO fix, 0.0103)
    
    % exchange function: power
    ep_power(1) = 16+rand;    %       * ep(1) = required daylight power, W                      (ALICINO transfer function)
    ep_power(2) = 16+rand;    %       * ep(2) = required eclipse power, W
    
    
    
    %----------------------------------------------------------------------
    tic
    for i=1:1000
        [M_power, P_power, info_power] = space_power_5subsystems(x_power,ep_power);
    end
    toc
    %----------------------------------------------------------------------
    
elseif run_flag == 6
    %======================================================================
    % SATELLITE
    %======================================================================
    d = rand(1,12);
    u = rand(1,21);
    par.fix.time = 365;
    
    %----------------------------------------------------------------------
    profile on
    tic
    for i=1:1000
        output = IAC2018_only_mass2_for_belief(d, u, par) ;
    end
    toc
    profile viewer
    %----------------------------------------------------------------------
    
    
elseif run_flag == 7
    %======================================================================
    % RES
    %======================================================================
    d = rand(1,12);
    u = rand(1,21);
    
    TIME_vector = 0:365/150:365;
    N_images_vector = ones(1,151);
    R1_function = @(x)(R1fun(TIME_vector,N_images_vector,x));
    R2_function = @(x)(R2fun(TIME_vector,N_images_vector,x));
    
    
    TM = 365;
    
    lam0 = 1;
    mu0 = 2;
    [lam,mu] = cubesat_failure_rate(d,u, lam0, mu0);
    lam = lam/365;
    mu = mu/365;
    
    %----------------------------------------------------------------------
    tic
    for i=1:1000
        RES = resilienceSimpleFunction(TM,R2_function,R1_function,lam, mu);
    end
    toc
    %----------------------------------------------------------------------
    
end