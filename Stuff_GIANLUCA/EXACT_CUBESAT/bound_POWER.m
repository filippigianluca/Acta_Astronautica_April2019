function [lb_d, ub_d, lb_u, ub_u] = bound_POWER()


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




lb_d = [ 0, ...   %       * x(1) = type of solar cell
         3, ...   %       * x(2) = required bus voltage, V [0,100]
         0, ...   %       * x(3) = eps configuration (det, mppt) [0,0.5) = det; [0.5,1] = mppt    
       0.5, ...   %       * x(4) = cell packing efficiency (or assembly factor) [0,1]
       0.5, ...   %       * x(5) = depth of discharge (it is how much the battery has been discharged) [0,1]
         1, ...   %       * x(6) = number of loads (to determine LCLs in PCDU) [1,50]
       0.4, ...   %       * x(7) = harness mass factor (as a percentage of eps mass) (harness=cables and buses)
         2, ...   %       * x(8) = allowable voltage drop, %
];



ub_d = [ 1, ...   %       * x(1) = type of solar cell
         5, ...   %       * x(2) = required bus voltage, V [0,100]
         1, ...   %       * x(3) = eps configuration (det, mppt) [0,0.5) = det; [0.5,1] = mppt    
       0.8, ...   %       * x(4) = cell packing efficiency (or assembly factor) [0,1]
         1, ...   %       * x(5) = depth of discharge (it is how much the battery has been discharged) [0,1]
        10, ...   %       * x(6) = number of loads (to determine LCLs in PCDU) [1,50]
       0.5, ...   %       * x(7) = harness mass factor (as a percentage of eps mass) (harness=cables and buses)
       2.5, ...   %       * x(8) = allowable voltage drop, %
];



lb_u = [ 0, ...   %       * ep(1) = required daylight power, W
         0, ...   %       * ep(2) = required eclipse power, W
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


ub_u = [ 30, ...   %       * ep(1) = required daylight power, W
         30, ...   %       * ep(2) = required eclipse power, W
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
return