function [lb_d, ub_d, lb_u, ub_u, par] = bound_SPACECRAFT()

%% space_aocs: Attitude and orbit control system model
%
% [M,P,varargout] = space_aocs(x,ep)
%
%% Inputs:
% * x: Design parameters
%       * x(1) = pointing accuracy, deg
%       * x(2) = offset between centre of gravity and centre of pressure, m
%       * x(3) = area normal to velocity vector, m^2
%       * x(4) = area normal to sun line, m^2
%       * x(5) = reflectance factor
%       * x(6) = spacecraft residual dipole, A*m^2
%       * x(7) = drag coefficient
%       * x(8) = actuator type: [0 0.5) = magtorq, [0.5 1] = thruster
%       * x(9) = specific impulse, s
%       * x(10) = burn time for momentum dumping, s
%       * x(11) = slew angle, deg
%       * x(12) = time for slew maneuvre, s
%       * x(13) = initial tumbling rate, rad/s
%       * x(14) = spin angle, deg      
% * ep: Environmental parameters
%       * ep(1) = orbit period, s
%       * ep(2) = mean planet's magnetic field stregth, T
%       * ep(3) = mean incident solar radiation, W/m^2
%       * ep(4) = mean dynamic pressure, Pa
%       * ep(5) = mean gravitational field strength, deg
%       * ep(6) = angle between spacecraft z axis and nadir vector, deg
%       * ep(7) = spacecraft moment of inertia about x axis, kg*m^2
%       * ep(8) = spacecraft moment of inertia about y axis, kg*m^2
%       * ep(9) = spacecraft moment of inertia about z axis, kg*m^2
%       * ep(10) = thruster moment arm, m
%       * ep(11) = number of thrusters
%       * ep(12) = number of dumping maneuvres during S/C lifetime
%       * ep(13) = number of slew maneuvres during S/C lifetime

lb_d = [ 20, ...  %       * x(1) = pointing accuracy, deg 
         0, ...   %       * x(2) (ALICINO d_aocs_1)= offset between centre of gravity and centre of pressure, m
         0.034, ...   %   * x(3) (ALICINO) = area normal to velocity vector, m^2  [0.034 0.15]
         0.034, ...   %   * x(4) (ALICINO) = area normal to sun line, m^2         
         0.5, ... %       * x(5) (ALICINO) = reflectance factor                   [0.5 0.7]           
         0.0005,  ...  %  * x(6) (ALICINO) = spacecraft residual dipole, A*m^2    [0.0005 0.0015]
         2, ...   %       * x(7) (ALICINO) = drag coefficient                     [2 2.5]
         0.1, ...   %     * x(8) = actuator type: [0 0.5) = magtorq, [0.5 1] = thruster
         100, ...   %     * x(9) = specific impulse, s
         10, ...   %      * x(10) = burn time for momentum dumping, s
         10, ...  %       * x(11)(ALICINO) = slew angle, deg
         30, ...  %       * x(12)(ALICINO)  = time for slew maneuvre, s
         10, ...   %      * x(13) = initial tumbling rate, rad/s
         -100 ...    %    * x(14) = spin angle, deg 
];

     
     
ub_d = [60, ...   %        * x(1) = pointing accuracy, deg
         0.3, ... %        * x(2) = offset between centre of gravity and centre of pressure, m 
         0.15, ...   %     * x(3) (ALICINO u_aocs_2) = area normal to velocity vector, m^2
         0.15, ...   %     * x(4) (ALICINO) = area normal to sun line, m^2
         0.7, ...   %      * x(5) (ALICINO) = reflectance factor
         0.0015, ...   %   * x(6) (ALICINO) = spacecraft residual dipole, A*m^2
         2.5, ...  %       * x(7) (ALICINO) = drag coefficient 
         0.3, ...  %       * x(8) = actuator type: [0 0.5) = magtorq, [0.5 1] = thruster
         200, ...  %       * x(9) = specific impulse, s
         50, ...  %       * x(10) = burn time for momentum dumping, s
         60, ...   %       * x(11) (ALICINO) = slew angle, deg
         90, ...   %       * x(12) (ALICINO) = time for slew maneuvre, s
         20, ...   %       * x(13) = initial tumbling rate, rad/s
         100 ...    %       * x(14) = spin angle, deg 
];

lb_u = [4000, ...       %       * ep(1) (ALICINO) = orbit period, s
        3e-6, ...       %       * ep(2) (ALICINO) = mean planet's magnetic field stregth, T
        500, ...        %       * ep(3) (ALICINO) = mean incident solar radiation, W/m^2
        1e-9, ...       %       * ep(4) = mean dynamic pressure, Pa
        1e-7, ...       %       * ep(5) = mean gravitational field strength, deg
        0, ...          %       * ep(6) = angle between spacecraft z axis and nadir vector, deg
        0.04, ...       %       * ep(7) (ALICINO) = spacecraft moment of inertia about x axis, kg*m^2
        0.01, ...       %       * ep(8) (ALICINO) = spacecraft moment of inertia about y axis, kg*m^2
        0.01, ...       %       * ep(9) (ALICINO) = spacecraft moment of inertia about z axis, kg*m^2
        1, ...          %       * ep(10) = thruster moment arm, m
        1, ...          %       * ep(11) = number of thrusters
        1, ...          %       * ep(12) = number of dumping maneuvres during S/C lifetime
        1, ...          %       * ep(13) = number of slew maneuvres during S/C lifetime
];



ub_u = [7000, ...      %       * ep(1) (ALICINO) = orbit period, s
        3e-3, ...      %       * ep(2) (ALICINO) = mean planet's magnetic field stregth, T
        2500, ...      %       * ep(3) (ALICINO) = mean incident solar radiation, W/m^2
        1e-7, ...        %     * ep(4) = mean dynamic pressure, Pa
        1e-3, ...       %      * ep(5) = mean gravitational field strength, deg
        30, ...       %       * ep(6) = angle between spacecraft z axis and nadir vector, deg
        1, ...         %       * ep(7) (ALICINO) = spacecraft moment of inertia about x axis, kg*m^2
        1, ...         %       * ep(8) (ALICINO) = spacecraft moment of inertia about y axis, kg*m^2
        1, ...         %       * ep(9) (ALICINO) = spacecraft moment of inertia about z axis, kg*m^2
        2, ...         %       * ep(10) = thruster moment arm, m
        5, ...         %       * ep(11) = number of thrusters
        5, ...         %       * ep(12) = number of dumping maneuvres during S/C lifetime
        15, ...         %       * ep(13) = number of slew maneuvres during S/C lifetime
];



par.dim_d_aocs = length(lb_d);
par.dim_u_aocs = length(lb_u);





%% space_ttc: Communications system model
%
%   [M,P,info] = space_ttc(x,ep)
%
%% Inputs:
% * x: Design parameters
%       * x(1) = frequency, GHz
%       * x(2) = modulation,  
%                   PSK  = [0 1/8)
%                   BPSK = [1/8 2/8) 
%                   CFSK = [2/8 3/8) 
%                   BFSK = [3/8 4/8)
%                   FSK  = [4/8 5/8)
%                   DPSK = [5/8 6/8)
%                   QPSK = [6/8 7/8)
%                   NRZ  = [7/8 1]
%       * x(3) = antenna efficiency
%       * x(4) = antenna gain, dB
%       * x(5) = onboard loss, dB
%       * x(6) = other unmodelled losses, dB (polarization,implementation...)
%       * x(7) = mass of distribution network, kg
%       * x(8) = modulation index, rad
%       * x(9) = amplifier type, TWTA = [0 0.5), SSPA = [0.5 1]
%
% * ep: Environmental parameters
%       * ep(1) = Bit Error Rate
%       * ep(2) = data volume, bits
%       * ep(3) = ground station G/T, dB
%       * ep(4) = range, km
%       * ep(5) = elevation angle, deg
%       * ep(6) = pointing accuracy, deg
%       * ep(7) = ground station antenna diameter, m
%       * ep(8) = ground station access time, min
%       * ep(9) = link margin, dB



lb_d = [ lb_d, ...    
         7, ...   %       * x(1) = frequency, GHz
         0, ...   %       * x(2) = modulation,  
       0.4, ...   %       * x(3) = antenna efficiency
         1, ...   %       * x(4) = antenna gain, dB
        0.1, ...  %       * x(5) = onboard loss, dB
        0.5, ...  %       * x(6) = other unmodelled losses, dB (polarization,implementation...)
        0.1, ...  %       * x(7) = mass of distribution network, kg
        0.5, ...  %       * x(8) = modulation index, rad
          0, ...  %       * x(9) = amplifier type, TWTA = [0 0.5), SSPA = [0.5 1]
];
     
     
ub_d = [ub_d, ...    
         10, ...   %       * x(1) = frequency, GHz
          1, ...   %       * x(2) = modulation,  
        0.6, ...   %       * x(3) = antenna efficiency
          5, ...   %       * x(4) = antenna gain, dB
          1, ...   %       * x(5) = onboard loss, dB
          2, ...   %       * x(6) = other unmodelled losses, dB (polarization,implementation...)
        0.5, ...   %       * x(7) = mass of distribution network, kg
          1, ...   %       * x(8) = modulation index, rad
          0.4         , ...   %       * x(9) = amplifier type, TWTA = [0 0.5), SSPA = [0.5 1]
];

     
lb_u = [lb_u, ...    
        1e-6, ...       %       * ep(1) = Bit Error Rate
         1e5, ...       %       * ep(2) = data volume, bits
          20, ...       %       * ep(3) = ground station G/T, dB
         600, ...       %       * ep(4) = range, km
           5, ...       %       * ep(5) = elevation angle, deg
           5, ...       %       * ep(6) = pointing accuracy, deg
          10, ...       %       * ep(7) = ground station antenna diameter, m
           5, ...       %       * ep(8) = ground station access time, min
           5, ...       %       * ep(9) = link margin, dB
];

     
ub_u = [ub_u, ... 
        1e-4, ...       %       * ep(1) = Bit Error Rate
         1e7, ...       %       * ep(2) = data volume, bits
          40, ...       %       * ep(3) = ground station G/T, dB
         700, ...       %       * ep(4) = range, km
          15, ...       %       * ep(5) = elevation angle, deg
          15, ...       %       * ep(6) = pointing accuracy, deg
          20, ...       %       * ep(7) = ground station antenna diameter, m
          15, ...       %       * ep(8) = ground station access time, min
          15, ...       %       * ep(9) = link margin, dB
];




par.dim_d_ttc = length(lb_d);
par.dim_u_ttc = length(lb_u);





%% space power: Electrical power system model
%
% [M,P,info] = space_power(x,ep)
%
%% Inputs:
% * x: Design parameters:
%       * x(1) = type of solar cell
%               Si = [0 0.33], GaAs28 = [0.33 0.66], GaAs30 = [0.66 1] 
%       * x(2) = required bus voltage, V [0,100]
%       * x(3) = eps configuration (det, mppt) [0,0.5) = det; [0.5,1] = mppt    
%              eps: electric power system, der: direct energy transfer, mppt: maximum point power tracking
%              (it is how the solar cells are connected to the battery, det=directly; mppt=using an interface that optimises the performance)
%       * x(4) = cell packing efficiency (or assembly factor) [0,1]
%       * x(5) = depth of discharge (it is how much the battery has been discharged) [0,1]
%       * x(6) = number of loads (to determine LCLs in PCDU) [1,50]
%              PCDU(Power Control and Distribution Unit):controls the power flow from the solar generator to the battery and distributes power, on command from the onboard computer, to the instruments, the heaters and the propulsion system
%              LCL(Latching Current Limiters)
%       * x(7) = harness mass factor (as a percentage of eps mass) (harness=cables and buses)
%       * x(8) = allowable voltage drop, %
%
% * ep: Environmental parameters:
%       * ep(1) = required daylight power, W
%       * ep(2) = required eclipse power, W
%       * ep(3) = orbital daylight time, h
%       * ep(4) = orbital eclipse time, h
%       * ep(5) = distance from Sun, AU
%       * ep(6) = worst case angle of incidence, deg
%       * ep(7) = uncertainty on array temperature, C
%       * ep(8) = satellite lifetime, yrs
%       * ep(9) = bus regulation (0: unregulated, 1: regulated)
%       * ep(10) = uncertainty on power requirements, %
%       * ep(11) = total radiation fluence
%       * ep(12) = s/c type (1: cubesat, 0: big one)




lb_d = [lb_d, ...
          0, ...   %       * x(1) = type of solar cell
         3, ...   %       * x(2) = required bus voltage, V [0,100]
         0, ...   %       * x(3) = eps configuration (det, mppt) [0,0.5) = det; [0.5,1] = mppt    
       0.5, ...   %       * x(4) = cell packing efficiency (or assembly factor) [0,1]
       0.5, ...   %       * x(5) = depth of discharge (it is how much the battery has been discharged) [0,1]
         1, ...   %       * x(6) = number of loads (to determine LCLs in PCDU) [1,50]
       0.4, ...   %       * x(7) = harness mass factor (as a percentage of eps mass) (harness=cables and buses)
         2, ...   %       * x(8) = allowable voltage drop, %
];



ub_d = [ub_d, ...
          1, ...   %       * x(1) = type of solar cell
         5, ...   %       * x(2) = required bus voltage, V [0,100]
         1, ...   %       * x(3) = eps configuration (det, mppt) [0,0.5) = det; [0.5,1] = mppt    
       0.8, ...   %       * x(4) = cell packing efficiency (or assembly factor) [0,1]
         1, ...   %       * x(5) = depth of discharge (it is how much the battery has been discharged) [0,1]
        10, ...   %       * x(6) = number of loads (to determine LCLs in PCDU) [1,50]
       0.5, ...   %       * x(7) = harness mass factor (as a percentage of eps mass) (harness=cables and buses)
       2.5, ...   %       * x(8) = allowable voltage drop, %
];



lb_u = [lb_u, ...
         1, ...   %       * ep(3) = orbital daylight time, h
         1, ...   %       * ep(4) = orbital eclipse time, h
       0.8, ...   %       * ep(5) = distance from Sun, AU
         0, ...   %       * ep(6) = worst case angle of incidence, deg
         0.5, ...   %       * ep(7) = uncertainty on array temperature, C
       0.5, ...   %       * ep(8) = satellite lifetime, yrs
         0.7, ...   %       * ep(9) = bus regulation (0: unregulated, 1: regulated)
       0.3, ...   %       * ep(10) = uncertainty on power requirements, %
         1e11, ...   %       * ep(11) = total radiation fluence
         0, ...   %       * ep(12) = s/c type (1: cubesat, 0: big one)
];


ub_u = [ ub_u, ...
          2, ...   %       * ep(3) = orbital daylight time, h
          2, ...   %       * ep(4) = orbital eclipse time, h
        1.2, ...   %       * ep(5) = distance from Sun, AU
         50, ...   %       * ep(6) = worst case angle of incidence, deg
         30, ...   %       * ep(7) = uncertainty on array temperature, C
        1.5, ...   %       * ep(8) = satellite lifetime, yrs
          1, ...   %       * ep(9) = bus regulation (0: unregulated, 1: regulated)
        0.6, ...   %       * ep(10) = uncertainty on power requirements, %
          1e16, ...   %       * ep(11) = total radiation fluence
          0.3, ...   %       * ep(12) = s/c type (1: cubesat, 0: big one)
];


par.dim_d_power = length(lb_d);
par.dim_u_power = length(lb_u);


return