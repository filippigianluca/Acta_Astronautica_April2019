%==========================================================================
% D E S C R I P T I O N:
%
% This scrip will run LambTAN_v1 in which the Edelbaum feasibility check 
% has been implemented. 
%==========================================================================


%==========================================================================
%                                                                         %
%    A D D     R E Q U I R E D     L I B R A R I E S                      %
%                                                                         %
%=========================================================================

clc
clear 

% Add the SpaceArt toolbox Library
addpath(genpath('astro_tool'))

%==========================================================================
%                                                                         %
%    I N I T I A L I Z E     P R O B L E M      P A R A M E T E R S       %
%                                                                         %
%==========================================================================

% State the Availabe Celestial Bodies
global_config.celestialBodies = { Planets.EARTH,         ... % The Initial Departure Plane----------------------t         
                                  Asteroids.neo163693,   ... % Asteroid #01
                                  Asteroids.neo164294,   ... % Asteroid #02 
                                  Asteroids.neo1998DK36, ... % Asteroid #03
                                  Asteroids.neo2004JG6,  ... % Asteroid #04
                                  Asteroids.neo2005TG45, ... % Asteroid #05
                                  Asteroids.neo2006WE4,  ... % Asteroid #06
                                  Asteroids.neo2007EB26, ... % Asteroid #07
                                  Asteroids.neo2008EA32, ... % Asteroid #08
                                  Asteroids.neo2008UL90, ... % Asteroid #09
                                  Asteroids.neo2010XB11, ... % Asteroid #10
                                  Asteroids.neo2012VE46, ... % Asteroid #11
                                  Asteroids.neo2013JX28};    % Asteroid #12
                                                       
                              
% global_config.celestialBodies = { Planets.EARTH,         ... % The Initial Departure Plane----------------------t         
%                                   Asteroids.neo163693,   ... % Asteroid #01
%                                   Asteroids.neo164294,   ... % Asteroid #02 
%                                   Asteroids.neo1998DK36, ... % Asteroid #03
%                                   Asteroids.neo2012VE46, ... % Asteroid #05
%                                   Asteroids.neo2005TG45, ... 
%                                   Asteroids.neo2006WE4,  ... % Asteroid #06
%                                   Asteroids.neo2007EB26};    % Asteroid #12
                                  
                                                                                         
                              
% Maximum allowed departure deltaV for Asteroids [km/s]
dVmax_asteroids = 1.5;  % original was 1.5
                              
% Maximum Allowed depature deltaV for each celestial body [km/s]
global_config.deltaVDepartur_Limit = dVmax_asteroids * ones(1, length(global_config.celestialBodies));    

% Set the Maximum departure deltaV from departure body [km/s]
global_config.deltaVDepartur_Limit(1) = 3;

% Minimum ToF [day]
global_config.tof_min = 30; 
                
% Maximum ToF [day]
global_config.tof_max = 365; 

% Time Step [day]
global_config.time_step = 10;

% Initial Problem Epoch [MJD2K]
global_config.epoch_start = date2mjd2000([2020,1,1,0,0,0]);

% End Problem Epoch [MJD2K]
global_config.epoch_end   = date2mjd2000([2030,1,1,0,0,0]);

% Global Maximum ToF [day]
global_config.tof_global_max = global_config.epoch_end - global_config.epoch_start; 

% Minimum Low-Thrust Acceleration [m/s^2]  
global_config.low_thrust_acc = 10^-4; 

% Mx's Coefficient to check feasibility 
global_config.low_trust_dVmax_coef = 2.0; 
                
% The Minimum allowed perihelium for a trajectory [km]
global_config.perihelium_min = 46001200; % Perihelium from Mercury

% Set Debug Flag
global_config.show_debug_level1 = false;

% Get Sun Gravitational Parameter [km^3/s^2] 
mu = AstroConstants.Sun_Planetary_Const;                                      

% Name of the output
output_file = 'solutions_20150712.mat';

NOTE = 'This test case is to gather solution feasible for both cases: Max (Coeff = 2.0) and  Endelbaum';
NOTE = [NOTE '\nThis simulation has been set in order to evaluate how much change the found solutions for a MX Coef. of 2. agains the previous simulation with MXC = 3.'];


% Global Solutions
global solutions;
global func_calls;
global num_solutions;
solutions     = {};
func_calls    = 0;
num_solutions = 0;

%==========================================================================
%                                                                         %
%    S T A R T     T H E    P R O B L E M     S O L U T I O N S           %
%                                                                         %
%==========================================================================

tic

% Initialize the current sequence
curr_seq = {}; 

% Initialize the missing celestial bodies
missing_CelestialBodies = 2:length(global_config.celestialBodies); 

% Get the Current Departure Celestial Body
curr_departure_celbody = global_config.celestialBodies{1};

% Run the Recursive Search and Sequence creation (LambTAN)
solutions2 = LambTAN_v1( curr_departure_celbody,  ... % Current Departure Celestial Body
                         1,                       ... % Index of the Departure Celestial Body
                         curr_seq,                ... % Current Trajectory Arc
                         missing_CelestialBodies, ... % Missing Celestial Bodies Array 
                         global_config,           ... % Global Problem Setting
                         mu);                     ... % Gravity Parameter
% Get the elapsted time                                       
time_elapsed = toc;

% Convert elapsed time from seconds to HH MM SS
[hh mm ss] = seconds2HHMMSS(time_elapsed);

% Display the Elapsed time
disp('-----------------------------------------------------------------------');
disp(sprintf('Elapsed time: %d hours %d minutes %f seconds.', hh, mm, ss));
disp(sprintf('Total Func Calls: %d.', func_calls ));
disp('-----------------------------------------------------------------------');


% Save the solutions
%solutions_sorted = sortSolution(solutions, 0.5, 0.5);
%save(output_file, 'global_config', 'solutions', 'time_elapsed', 'solutions_sorted', 'NOTE');
