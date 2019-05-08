function [F] = Spacecraft(d, u, ~)  % dim(d)=10    dim(u)=16; 
%
% 15 u parameters, new order after christmass
%
% INPUT:
%
% u_1:  AOCS l [m]
% u_2:  AOCS A [m^2]
% u_3:  AOCS q []
% u_4:  AOCS m [mA*m^2]
% u_5:  AOCS C_D []
% u_6:  AOCS dI [%]

% u_7:  TTC eta_ant []
% u_8:  TTC G_t [dB]
% u_9:  TTC L_t [dB]
% u_10: TTC L_other [dB]
% u_11: TTC M_rfdn  [kg]

% u_12: EPS D_cell []
% u_13: EPS eta_a []
% u_14: EPS rho_sa [kg/m^2]
% u_15: EPS dP [%]
% u_16: EPS T_max [C]





global nfevalglobal;
        nfevalglobal = nfevalglobal + 1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  before  order  vector
%
% u(1)= u_sensitivity(8); %u_sensitivity(1);
% u(2)= u_sensitivity(9); %u_sensitivity(2);
% u(3)= u_sensitivity(10); %u_sensitivity(3);
% u(4)= u_sensitivity(12);
% u(5)= u_sensitivity(11); %u_sensitivity(4);
% u(6)= u_sensitivity(13);
% u(7)= u_sensitivity(1); %u_sensitivity(5);
% u(8)= u_sensitivity(14);
% u(9)= u_sensitivity(15);
% u(10)=u_sensitivity(16);
% u(11)=u_sensitivity(2);%u_sensitivity(6);
% u(12)=u_sensitivity(3);%u_sensitivity(7);
% u(13)=u_sensitivity(4);%u_sensitivity(8);
% u(14)=u_sensitivity(5);%u_sensitivity(9);
% u(15)=u_sensitivity(6);%u_sensitivity(10);
% u(16)=u_sensitivity(7);%u_sensitivity(11);
% 
% 
% 
% 
% % u(1)=u_sensitivity(1);
% % u(2)=u_sensitivity(2);
% % u(3)=u_sensitivity(3);
% % u(4)=u_sensitivity(12);
% % u(5)=u_sensitivity(4);
% % u(6)=u_sensitivity(13);
% % u(7)=u_sensitivity(5);
% % u(8)=u_sensitivity(14);
% % u(9)=u_sensitivity(15);
% % u(10)=u_sensitivity(16);
% % u(11)=u_sensitivity(6);
% % u(12)=u_sensitivity(7);
% % u(13)=u_sensitivity(8);
% % u(14)=u_sensitivity(9);
% % u(15)=u_sensitivity(10);
% % u(16)=u_sensitivity(11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% AOCS
%%%%%%%%%%%%%%%%%%%
x_aocs(1)  = 5;                      % FIXED: pointing accuracy, deg
ep_aocs(1) = 5850;                   % FIXED: orbit period, s
x_aocs(8)  = 0.1;                    % FIXED: actuator type: [0 0.5) = magtorq, [0.5 1] = thruster
x_aocs(9)  = 1;                      % FIXED: specific impulse, s is about thruster
x_aocs(10) = 1;                      % FIXED: burn time for momentum dumping, s
x_aocs(13) = 0;                      % FIXED: initial tumbling rate, rad/s
x_aocs(14) = 0;                      % FIXED: spin angle, deg 
ep_aocs(2) = 3.1e-5;                 % FIXED: mean planet's magnetic field stregth, T
ep_aocs(3) = 1420;                   % FIXED: mean incident solar radiation, W/m^2
ep_aocs(4) = 0.5*1e-8*(7.54^2);      % FIXED: mean dynamic pressure, Pa      (1/2 rho v^2)
ep_aocs(5) = 3*398600/(6378 +640)^3; % FIXED: mean gravitational field strength, deg
ep_aocs(6) = 5;                      % FIXED: angle between spacecraft z axis and nadir vector, deg
ep_aocs(10) = 0;                     % FIXED: thruster moment arm, m
ep_aocs(11) = 0;                     % FIXED: number of thrusters
ep_aocs(12) = 0;                     % FIXED: number of dumping maneuvres during S/C lifetime
ep_aocs(13) = 0;                     % FIXED: number of slew maneuvres during S/C lifetime



x_aocs(2) = u(1);                    % offset between centre of gravity and centre of pressure, m
x_aocs(3) = u(2);                    % area normal to velocity vector, m^2
x_aocs(4) = u(2);                    % area normal to sun line, m^2
x_aocs(5) = u(3);                    % reflectance factor
x_aocs(6) = u(4)/1000;               % spacecraft residual dipole, A*m^2   
x_aocs(7) = u(5);                    % drag coefficient

ep_aocs(7) = 0.0417*(1+u(6)/100);    % spacecraft moment of inertia about x axis, kg*m^2
ep_aocs(8) = 0.1083*(1+u(6)/100);    % spacecraft moment of inertia about y axis, kg*m^2
ep_aocs(9) = 0.1417*(1+u(6)/100);    % spacecraft moment of inertia about z axis, kg*m^2




x_aocs(11)= d(1);                    % slew angle, deg
x_aocs(12)= d(2);                    % time for slew maneuvre, s


         
[M_aocs, P_aocs, ~] = space_aocs(x_aocs,ep_aocs);
%%%%%%%%%%%%%%%%%%%%%%%%


%%TTC
%%%%%%%%%%%%%%%%%%%%%%%%
x_ttc(1) = d(3);              % frequency, GHz
x_ttc(2) = d(4);              % modulation
x_ttc(3) = u(7);              % antenna efficiency
x_ttc(4) = u(8);              % antenna gain, dB
x_ttc(5) = u(9);              % onboard loss, dB
x_ttc(6) = u(10);             % other unmodelled losses, dB (polarization,implementation...)
x_ttc(7) = u(11);             % mass of distribution network, kg
x_ttc(8) = 1;                 % FIXED: modulation index, rad (beta)
x_ttc(9) = d(5);              % amplifier type, TWTA = [0 0.5), SSPA = [0.5 1]


ep_ttc(1) = 1e-5;             % FIXED: Bit Error Rate             
ep_ttc(2) = 1e6;              % FIXED: data volume, bits
ep_ttc(3) = 30;               % FIXED: ground station G/T, dB/K 
ep_ttc(4) = 640;              % FIXED: range, km
ep_ttc(5) = 10;               % FIXED: elevation angle, deg
ep_ttc(6) = ep_aocs(6);       % FIXED: pointing accuracy, deg
c = 299792458;
lambda = c/(1e9*x_ttc(1));
ep_ttc(7) = lambda/pi*sqrt(10^(6)/0.55);
ep_ttc(8) = 600/60;           % FIXED
ep_ttc(9) = 0;                % FIXED

[M_ttc, P_ttc, ~] = space_ttc(x_ttc, ep_ttc);
%%%%%%%%%%%%%%%%%%%%%%%%%



%%POWER
%%%%%%%%%%%%%%%%%%%%%%%%%
CL = 14250;                              % FIXED:   number of charge/discharge
x_power(1)  = 0.7;                       % FIXED:   type of solar cell  (Si = [0 0.33], GaAs28 = [0.33 0.66], GaAs30 = [0.66 1])  
ep_power(3) = 1.615;                     % FIXED:   orbital daylight time, h
ep_power(4) = 1;% 0.0103;%                   % FIXED:   orbital eclipse time, h
ep_power(5) = 1;                         % FIXED:   distance from Sun, AU
ep_power(6) = 5;                         % FIXED:   worst case angle of incidence, deg
ep_power(7) = u(16);                    % MODIFIED:    margin T_max   uncertainty on array temperature, C
ep_power(8) = 1;                         % FIXED:   satellite lifetime, yrs
ep_power(9) = 1;                         % FIXED:   bus regulation (0: unregulated, 1: regulated)
ep_power(11)= 1;                         % FIXED:   total radiation fluence
ep_power(12)= 1;                         % FIXED:   s/c type (1: cubesat, 0: big one)



x_power(2) = d(8);                       % required bus voltage, V [0,100]           
x_power(3) = d(10);                      % eps configuration (det, mppt) [0,0.5) = det; [0.5,1] = mppt           
x_power(4) = u(13);  %0.8107;%                    % cell packing efficiency (or assembly factor) [0,1]
x_power(5) = -36.76*log(CL/207800)/100;  % DOD: depth of discharge [0,1] (Alicino e Vasile paper)
x_power(6) = 0;%??             
x_power(7) = 0;%??          
x_power(8) = d(9);                       % allowable voltage drop (%)   (x_power(2)*d(9);)    

ep_power(1) = 16+P_ttc+10*P_aocs;           % required daylight power, W                     
ep_power(2) = 16+P_ttc+10*P_aocs;           % required eclipse power, W
ep_power(10) = u(15);                    % uncertainty on power requirements, % 



x_power(9) = u(14);                      % MODIFIED:    rho_sa, solar array density
ep_power(13) = d(6);                     % MODIFIED:    eta_cell, cell efficiency
ep_power(14) = d(7);                     % MODIFIED:    E_cell, Wh
ep_power(15) = u(12);                    % MODIFIED:    degradation
% ep_power(16) = u(16);                    % MODIFIED:    T_max


[M_power, ~, ~] = space_power(x_power,ep_power);
%%%%%%%%%%%%%%%%%%%%%%%%%


F = M_aocs + M_ttc + M_power;
end