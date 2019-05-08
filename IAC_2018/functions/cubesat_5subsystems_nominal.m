function [d, u] = cubesat_5subsystems_nominal()
% nominal valoues for the cubesat function with 5 subsystems: AOCS, TTC,
% POWER, PAYLOAD and OBDH.






%% PAYLOAD
%---design-----------------------------------------------------------------
d(10) = 10000;  % [0 10000]  TIME, fixed (sec)
d(11) = 4;      % [1 5]      bit_depth 
d(12) = 0.4;    % [0 1]      camera_catalogs 0.4

%---epistemic--------------------------------------------------------------
u(17) = 800;    % altitude (km)
u(18) = 0.2;    % elevation angle (deg)
u(19) = 0.2;    % delta inclination orbit (%)




%% AOCS 
%---design-----------------------------------------------------------------
% x_aocs(1) = 30;         %%%       * x(1) = pointing accuracy, deg                                                (fix, 30)   
u(1) = 0.01;       %       * x(2) = offset between centre of gravity and centre of pressure, m   (ALICINO uncertain, [0.005 0.01][0.01 0.02])
u(2) = 0.0885;     %       * x(3) = area normal to velocity vector, m^2                          (ALICINO uncertain, [0.034 0.0885][0.0885 0.15])
% x_aocs(4) = 0.0885;     %       * x(4) = area normal to sun line, m^2                                 (ALICINO uncertain, [0.034 0.0885][0.0885 0.15])
u(3) = 0.6;        %       * x(5) = reflectance factor                                           (ALICINO uncertain, [0.5 0.6][0.6 0.7])
u(4) = 0.001;      %       * x(6) = spacecraft residual dipole, A*m^2                            (ALICINO uncertain, [0.0005 0.001][0.001 0.0015])
u(5) = 2.2;        %       * x(7) = drag coefficient                                             (ALICINO uncertain, [2 2.2][2.2 2.5])
% x_aocs(8) = 0.1;        %%%       * x(8) = actuator type: [0 0.5) = magtorq, [0.5 1] = thruster                    (magnetotorque or thruster, fix, 0.1)          
% x_aocs(9) =  200;       %%%       * x(9) = specific impulse, s                                                   (only for thruster, fix 200) 
% x_aocs(10) =  50;       %%%       * x(10) = burn time for momentum dumping, s                                    (only for thruster, fix 50) 
d(1) =  10;       %       * x(11) = slew angle, deg                                             (ALICINO design, [10 60])
d(2) =  90;       %       * x(12) = time for slew maneuvre, s                                   (ALICINO design, [30 90])
% x_aocs(13) =  0;        %%%       * x(13) = initial tumbling rate, rad/s                                         (fix, 10)
% x_aocs(14) =  0;        %%%       * x(14) = spin angle, deg                                                      (fix, 100)

%---epistemic--------------------------------------------------------------
% ep_aocs(1) =  5850;     %%%       * ep(1) = orbit period, s                                             (ALICINO fixed, 5850)           
% ep_aocs(2) =  3.1e-5;   %%%       * ep(2) = mean planet's magnetic field stregth, T                     (ALICINO fixed, 3.1e-5)
% ep_aocs(3) =  1420;     %%%       * ep(3) = mean incident solar radiation, W/m^2                        (ALICINO fixed, 1420)
% ep_aocs(4) =  1/2*1e-5*7.54^2;%%%       * ep(4) = mean dynamic pressure, Pa                                      (1/2 * rho * v^2 = 1/2 * 1e-5  * 7.54^2)
% ep_aocs(5) = 1e-5;      %%%       * ep(5) = mean gravitational field strength, deg                               (fixed)
% ep_aocs(6) =  5;        %%%       * ep(6) = angle between spacecraft z axis and nadir vector, deg       (ALICINO fixed, 5)
% ep_aocs(7) =  0.0417;   %%%       * ep(7) = spacecraft moment of inertia about x axis, kg*m^2           (ALICINO fixed, 0.0417)
% ep_aocs(8) =  0.1083;   %%%       * ep(8) = spacecraft moment of inertia about y axis, kg*m^2           (ALICINO fixed, 0.1083)
% ep_aocs(9) =  0.1417;   %%%       * ep(9) = spacecraft moment of inertia about z axis, kg*m^2           (ALICINO fixed, 0.1417)
% ep_aocs(10) =  1;       %%%       * ep(10) = thruster moment arm, m                                              (only for thruster, fix 1) 
% ep_aocs(11) =  1;       %%%       * ep(11) = number of thrusters                                                 (only for thruster, fix 1) 
% ep_aocs(12) =  1;       %%%       * ep(12) = number of dumping maneuvres during S/C lifetime                     (fix 1) 
% ep_aocs(13) =  1;       %%%       * ep(13) = number of slew maneuvres during S/C lifetime                        (fix 1) 
u(6) = 0;
%---FUNCTION---------------------------------------------------------------
% [M_aocs, P_aocs, ~] = space_aocs(x_aocs, ep_aocs);




% nominal values for TTC 
%---design-----------------------------------------------------------------
d(3) = 8.253;       %       * x(1) = frequency, GHz                           (ALICINO design       [7 10])
d(4) = 1;           %       * x(2) = modulation,                              (ALICINO design       [0 1])
u(7) = 0.8;         %       * x(3) = antenna efficiency                       (ALICINO uncertainty  [0.6 0.8] [0.8 0.9])
u(8) = 3;           %       * x(4) = antenna gain, dB                         (ALICINO uncertainty  [1 3][3 5] trasmit AG)
u(9) = 0.5;         %       * x(5) = onboard loss, dB                         (ALICINO uncertainty  [0.1 0.5][0.5 1])
u(10) = 2;           %       * x(6) = other unmodelled losses, dB              (ALICINO uncertainty  [0.5 1.5][1.5 2])
u(11) = 0.3;         %       * x(7) = mass of distribution network, kg         (ALICINO uncertainty  [0.1 0.3][0.2 0.5])
% x_ttc(8) = 0.5;         %%%       * x(8) = modulation index, rad
d(5) = 0;           %       * x(9) = amplifier type,                          (ALICINO design       [0 1])
%---epistemic--------------------------------------------------------------
% ep_ttc(1) = 1e-5;       %%%       * ep(1) = Bit Error Rate                          (ALICINO fix, 1e-5)
% ep_ttc(2) = 1e6;        %%%       * ep(2) = data volume, bits                       (ALICINO fix, 1e6)
% ep_ttc(3) = 30;         %%%       * ep(3) = ground station G/T, dB                  (ALICINO fix, 30)
% ep_ttc(4) = 640;        %%%       * ep(4) = range, km                               (ALICINO fix, 640)
% ep_ttc(5) = 10;         %%%       * ep(5) = elevation angle, deg                    (ALICINO fix, 10)
% ep_ttc(6) = 5;          %%%       * ep(6) = pointing accuracy, deg                  (ALICINO fix, 5)
% c = 299792458;
% ep_ttc(7) = (c/(1e9*x_ttc(1)))/pi*(10e6/0.55)^0.5;     %%%       * ep(7) = ground station antenna diameter, m      (lambda/pi * (10e6/0.55)^0.5; lambda = c/(1e9*f))
% ep_ttc(8) = 10;         %%%       * ep(8) = ground station access time, min         (ALICINO fix, 10)
% ep_ttc(9) = 0;          %%%       * ep(9) = link margin, dB                         (margin on the P, i can fix to 0)
%---FUNCTION---------------------------------------------------------------
% [M_ttc, P_ttc, ~] = space_ttc(x_ttc,ep_ttc);



% nominal values for POWER 
%--------------------------------------------------------------------------
d(6) =  0;         %       * x(1) = type of solar cell                                                                                       (design, [0 1])
d(7) =  5;         %       * x(2) = required bus voltage, V [0,100]                                                                  (ALICINO design, [3 5])                                                                                             
d(8) =  1;         %       * x(3) = eps configuration (det, mppt) [0,0.5) = det; [0.5,1] = mppt                                      (ALICINO design, [0 1])
u(12) =  0.85;      %       * x(4) = cell packing efficiency (or assembly factor) [0,1]                                               (ALICINO uncertain, [0.8 0.85][0.85 0.9])   
% x_power(5) =   -36.76*log(14250/207800)/100;  %%%       * x(5) = depth of discharge (it is how much the battery has been discharged) [0,1]        (ALICINO fix, DOD = -36.76 ln(CL/207800), CL 14250)
% x_power(6) =   1;        %%%       * x(6) = number of loads (to determine LCLs in PCDU) [1,50]
% x_power(7) =   10;       %%%       * x(7) = harness mass factor (as a percentage of eps mass) (harness=cables and buses)
d(9) =   1;        %       * x(8) = allowable voltage drop,                                                                          (ALICINO design, [1 3])              

% ep_power(1) = 16+P_aocs+P_ttc; %       * ep(1) = required daylight power, W                      (ALICINO transfer function)
% ep_power(2) = 16+P_aocs+P_ttc; %       * ep(2) = required eclipse power, W                       (ALICINO transfer function)
% ep_power(3) = 1.615;     %%%       * ep(3) = orbital daylight time, h                              (ALICINO fix, 1.615)
% ep_power(4) = 0.0103;    %%%       * ep(4) = orbital eclipse time, h                               (ALICINO fix, 0.0103)
% ep_power(5) = 1;         %%%       * ep(5) = distance from Sun, AU                                 (ALICINO fix, 1)
% ep_power(6) = 5;         %%%       * ep(6) = worst case angle of incidence, deg                    (ALICINO fix, 5)
u(13) = 10;        %       * ep(7) = uncertainty on array temperature, C                   (ALICINO uncertain, [0 10][10 15])
% ep_power(8) = 1;         %%%       * ep(8) = satellite lifetime, yrs                               (ALICINO fix, 1)
% ep_power(9) = 1;         %%%       * ep(9) = bus regulation (0: unregulated, 1: regulated)         (ALICINO fix, 1)
u(14) = 10;       %       * ep(10) = uncertainty on power requirements, %                 (ALICINO uncertain, [0 10][10 20])
% ep_power(11) = 1e14;     %%%       * ep(11) = total radiation fluence                                      (fixed)
% ep_power(12) = 1;        %%%      * ep(12) = s/c type (1: cubesat, 0: big one)                     (ALICINO fix, 1)
%-- add to the function "space_power_Alicino_paper":
u(15) = 0;       %        uncertain         Delta rho_sa [% kg/m^2] density solar array, [0 15][15 30]
u(16) = 0;       %        uncertain         Delta D_cell [%] degradation solar array,    [0 50][50 100]
%--------------------------------------------------------------------------
% [M_power, P_power, ~] = space_power_Alicino_paper(x_power,ep_power);





%--------------------------------------------------------------------------
% sum of the masses of Attitude and Orbit Control, Tellecomunication and
% Power subsystems
%--------------------------------------------------------------------------
% M = M_aocs + M_ttc + M_power;

return