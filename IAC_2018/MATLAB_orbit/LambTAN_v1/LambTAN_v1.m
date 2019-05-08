function [ solutions2 ] = LambTAN_v1( curr_dep,                ... % Current Departure Celestial Body
                                      ind_curr_dep,            ... % index of the Departure Celestial Body
                                      prev_seq,                ... % Current Trajectory Arc
                                      missing_CelestialBodies, ... % Missing Celestial Bodies Array  
                                      global_config,           ... % Global Problem Setting
                                      mu)                      ... % Gravity Parameter
%--------------------------------------------------------------------------
% LambTAN (Lambert To Target at Nodel points)
%
% INPUTS:
%    * curr_dep                  Current Departure Celestial Body
%    * curr_seq                  Current Trajectory Arc. ( curr_seq = -1 if
%                                not prev arc has been created)
%    * missing_CelestialBodies   Missing Celestial Bodies Array
%    * global_config             Global Problem Setting
%         - celestialBodies         All the celestial bodies                       
%    * solutions                 Global Solution Array [OUTPUT]
%    * mu                        Gravity Parameter       
% 
% OUTPUTS:
%    * solutions               Status Flag
%
% CHANGES:
%  - [2014/05/06] This is the last version of the Algorithm used for the 
%  Workshop on May. 
%  - [2014/07/09] The Max's low thrust feasibility method has been changed by
%  the Edelbaum method. (Please refer to Orbital Mechanics by Chovotov (pag
%  344).
%  - [2016/07/05] Funtion Evaluation counter has been added.
%  - [2016/07/25] The function is renamed from SearchAndCreateNewArc_v4 to 
%  LambTAN_v1
%--------------------------------------------------------------------------                                              

% Global Store Solution Array
global solutions;

% Global Number of Function Calls
global func_calls;

% Number of found Solutions
global num_solutions; 

% Initialize the output solutions
solutions2 = {};

% Get the Current Departure Celestial Body
curr_departure_celbody = curr_dep; 
                                     
% Get the Total Time of Flight so far
% DELETEME -> if (~isstruct(curr_seq) && curr_seq == -1)
if isempty(prev_seq)
    total_tof = 0;
else
    total_tof = prev_seq{end}.total_tof;
end

% Loop through out all the missing celestial bodies
for ind_curr_missing_celbody = 1 : length(missing_CelestialBodies)

    % Get the index of the current target celestial body
    ind_curr_target_celbody = missing_CelestialBodies(ind_curr_missing_celbody);
        
    % Get the Current Target Celestial Body
    curr_target_celbody = global_config.celestialBodies{ind_curr_target_celbody}; 
             
    % Get the current Missing Celestial Bodies Array without the current
    % target celestial body.
    curr_missing_celestialbodies = missing_CelestialBodies(missing_CelestialBodies ~= ind_curr_target_celbody);
                  
    % Compute the Nodal Points [deg]                
    [Mnode1, Mnode2, error_status] = computeNodalPoints_M0(curr_departure_celbody, ... % The Departure Celestial Body
                                                           curr_target_celbody,    ... % The Arrival Celestial Body
                                                           mu);                    ... % Gravitational Parameter                                                         
    % Sanity Check
    if error_status == 1
        error('An error was found during the Nodal Point computation');        
    end                

    % Show Debug Info: LEVEL 1
    if global_config.show_debug_level1
        if isa(curr_departure_celbody, 'double')
            disp([sprintf('Earth vs %s ----------------------', curr_target_celbody.name)])    
        else 
            disp(sprintf('%-10s vs %s ----------------------', curr_departure_celbody.name, curr_target_celbody.name))    
        end
        disp([sprintf('   --> Mnode #1: %f\tMnode #2: %f', Mnode1, Mnode2)]);
    end
        
    % Set the Start of the Local Time Window [MDJ2000]         
    if isempty(prev_seq)
        curr_local_windows_epoch_start = global_config.epoch_start;
    else
        curr_local_windows_epoch_start = prev_seq{end}.arrival_epoch;
    end
      
    % Set the End of the Local Time Window [MDJ2000]
    curr_local_windows_epoch_end = global_config.epoch_end; 
    
    %======================================================================
    % W A R N I N G ! ! ! !
    %======================================================================     
    % It is assumed that all the target celestial bodies will be
    % CelestialObject instances and not double/integer (Planet Id). Only
    % the initial departure celestial body is considered to be a
    % double/integer (Planet Id). This is because we are using the current  
    % SpaceArt Toolbox to get the ephemerides of the planets. So the next 
    % curr_target_celbody is assummed to be a CelestialBody object and it
    % is expected to have the getKepElements() method. 
    %======================================================================
    
    % Compute the Epoch at the first Nodal point  
    [curr_pass_epochs_Mn1, error_status] = computeNodalPassesEpochs( curr_target_celbody.getKeplerianElements(), ... % Keplerian Orbital Elements
                                                                     Mnode1,                               ... % Mean Anomaly [deg]
                                                                     curr_local_windows_epoch_start,       ... % Start Local Windows Epoch [MDJ2000]
                                                                     curr_local_windows_epoch_end,         ... % End   Local Windows Epoch [MDJ2000]
                                                                     mu);                                  ... % Graviational Parameter

    % Sanity Check
    if error_status == 1
        error('An error was found during the First Nodal Point Epoch passes computation');        
    end  
        
    % Compute the Epoch at Nodal Time ( keporb, M0, nrev, mu)    
    [curr_pass_epochs_Mn2, error_status] = computeNodalPassesEpochs( curr_target_celbody.getKeplerianElements(), ... % Keplerian Orbital Elements
                                                                     Mnode2,                               ... % Mean Anomaly [deg]
                                                                     curr_local_windows_epoch_start,       ... % Start Local Windows Epoch [MDJ2000]
                                                                     curr_local_windows_epoch_end,         ... % End   Local Windows Epoch [MDJ2000]
                                                                     mu);                                  ... % Graviational Parameter

    % Sanity Check
    if error_status == 1
        error('An error was found during the Second Nodal Point Epoch passes computation');        
    end 
           
    % Concatenate the Epochs [MDJ2000]
    curr_pass_epochs = [curr_pass_epochs_Mn1 curr_pass_epochs_Mn2];    
        
    % Go Through out all the pass Epochs of the Node Mnode1     
    for ind_epoch = 1 : length(curr_pass_epochs)
    
        % Get the Current Pass Epoch
        curr_pass_epoch = curr_pass_epochs(ind_epoch);
                        
        % Compute the current Local Start epoch
        curr_epoch_start = curr_pass_epoch - global_config.tof_max;

        % Sanity Check: We are inside the Local Time Domain
        if curr_epoch_start < curr_local_windows_epoch_start
            curr_epoch_start = curr_local_windows_epoch_start;
        end        
        
        % Compute the current Local End epoch
        curr_epoch_end = curr_pass_epoch - global_config.tof_min;

        % Sanity Check: We are outside of the Local Time Domain, we jump
        if curr_epoch_end < curr_local_windows_epoch_start
            continue;
        end
        
        % Compute the position [km] and velocity [km/s] vectors at Arrival (Nodal Point)                
        [arrival_r, arrival_v] = StardustTool.CartesianElementsAt( curr_target_celbody, curr_pass_epoch);    
       
        % Loop Through out all the posible trajectories arcs
        for curr_departure_epoch = curr_epoch_start : global_config.time_step : curr_epoch_end
                              
            % Compute the new ToF for the trajectory [day]
            curr_tof = (curr_pass_epoch - curr_departure_epoch);
            
            % Sanity Check: Total ToF inside the ToF limit
            if (total_tof + curr_tof) > ( global_config.tof_global_max)
                continue;
            end
            
            % Get Position and Velocity of the DEPARTURE Celestial Body for
            % the current departure epoch
            [departure_r, departure_v] = StardustTool.CartesianElementsAt( curr_departure_celbody, curr_departure_epoch );
            
            % Compute the Lambert
            [~,~,~,err,vel_initial,vel_final,~,~] = lambertMR(departure_r,      ... % Initial Vector Position
                                                              arrival_r,        ... % Final position vector
                                                              curr_tof * 86400, ... % Time of flight [seconds]
                                                              mu,               ... % Planetary constant of the planet (mu = mass * G) [L^3/T^2]
                                                              0,                ... % Logical variable defining whether transfer is
                                                                                ... %   0: direct transfer from R1 to R2 (counterclockwise)
                                                                                ... %   1: retrograde transfer from R1 to R2 (clockwise)
                                                              0,                ... % Number of revolutions (Nrev):
                                                                                ... %   if Nrev = 0 ZERO-REVOLUTION transfer is calculated
                                                                                ... %   if Nrev > 0 two transfers are possible. Ncase should be
                                                                                ... %   defined to select one of the two.
                                                              0,                ... % Logical variable defining the small-a or large-a option in
                                                                                ... % case of Nrev>0:
                                                                                ... %   0: small-a option
                                                                                ... %   1: large-a option
                                                              2);                   % LambertMR options:
                                                                                    %   optionsLMR(1) = display options:
                                                                                    %     - 0: no display
                                                                                    %     - 1: warnings are displayed only when the algorithm does not converge
                                                                                    %     - 2: full warnings displayed
           
            % Increment the Function Calls counter
            func_calls = func_calls + 1;                                      
            
            % Sanity Check
            if err ~= 0
                
                % Warn the user but not break the execution
                warning('Some error was found during the LamberMR computation');                                
                
                % Discart this arc
                continue;                
            end
                        
            % Compute the Total DeltaV
            dv1(1) = vel_initial(1) - departure_v(1);
            dv1(2) = vel_initial(2) - departure_v(2);
            dv1(3) = vel_initial(3) - departure_v(3);
        
            dv2(1) = arrival_v(1) - vel_final(1);
            dv2(2) = arrival_v(2) - vel_final(2);
            dv2(3) = arrival_v(3) - vel_final(3);
        
            deltaV_Departure =  abs(norm(dv1));
            deltaV_Arrival   =  abs(norm(dv2));
            deltaV_Total     = (deltaV_Departure + deltaV_Arrival);
                
            % Sanity Check: Maximum deltaV @ Departure
            if deltaV_Departure > global_config.deltaVDepartur_Limit(ind_curr_dep);                
                continue;
            end
                    
            % Edelbaum's feasibility Check 
            deltaV_edelbaum = sqrt( norm(departure_v)^2 + norm(vel_final)^2 - 2*norm(departure_v)*norm(vel_final)); % [km/S]
            deltaV_max      = ( global_config.low_trust_dVmax_coef * deltaV_Departure); % [km/S]
            
            if deltaV_edelbaum > deltaV_max
                deltaV_limit   = deltaV_edelbaum;
            else 
                deltaV_limit   = deltaV_max;
            end
            
            % Compute maximum dv provided by the Low-Thrust engine [m/s]
            dvlowthrust = global_config.low_thrust_acc * curr_tof * 86400;
        
            % Sanity Check: Low-Thrust Constrain            
            if dvlowthrust < ( deltaV_limit * 1000) % [m/s]                      
                continue;
            end                
            
            % ==========================================================================================                                                                        
        
            % Compute the Keplerian elements of the transfer orbit at the departure position (a in km and angles in rad)               
            %[MARILENA]: comment it to test if any changes
%           kep_transfer_orbit = cart2kep([arrival_r, vel_final],mu);
            kep_transfer_orbit = cart2kep([departure_r, vel_initial],mu);
            
            % Eccentricity of transfer orbit
            ecc = kep_transfer_orbit(2);
            
            if ecc > 1
                error('Transfer orbit eccentricity greater than 1')
            end

            % True anomaly at departure position over transfer orbit [rad]
            theta = kep_transfer_orbit(6);

            % Compute the eccentric anomaly at the departure position [rad]
            cos_E = ( cos(theta) + ecc ) / (1 + ecc * cos(theta) );
            sin_E = ( sin(theta) * sqrt(1 - ecc^2) ) / (1 + ecc * cos(theta) );

            E = atan2(sin_E, cos_E);
            E = mod(E, 2*pi);
            
            % Compute the mean anomaly at the departure position [rad]
            M = E - ecc*sin(E); 
            
            % Mean anomaly at the departure position [deg]
            M = M * 180/pi;
            
            % Creat object to use as 1st input in SearchAndCreateNewArc
            % (angles in degrees and distance in AU!)
            au2km = AstroConstants.Astronomical_Unit;
            curr_departure_orbit = CelestialBody('transfer_orbit',             ... % Name of the CelestialBody (This case current trajectory
                                                 kep_transfer_orbit(1)/au2km,  ... % Semimajor axis [AU]
                                                 kep_transfer_orbit(2),        ... % Eccentricity
                                                 kep_transfer_orbit(3)*180/pi, ... % Inclination [deg]
                                                 kep_transfer_orbit(4)*180/pi, ... % Asc. Node/Raan [deg]
                                                 kep_transfer_orbit(5)*180/pi, ... % Arg. Perigee [deg]
                                                 M,                            ... % Mean anomoly, M at time given t0 [deg]
                                                 curr_departure_epoch ); % Time at which Mo is given [MJD2000]  
                                             
            % ==========================================================================================   
                               
            % Compute the Perihelium
            p = kep_transfer_orbit(1) * (1 - kep_transfer_orbit(2));
                        
            % Sanity check: too close to the sun
            if p < global_config.perihelium_min             
                continue; 
            end
                        
            % Sanity Check: First Arc            
            if isempty(prev_seq)     
                
                %disp('Creating the first ARC _____________________________________')
                
                % Create the First Arc
                arc = TrajectoryArc;
                arc.id               =  1;                       % Id for the arc as well to control the size of the sequence
                arc.deltaV_departure = deltaV_Departure;         % The deltaV @ Departure [km/s]
                arc.deltaV_arrival   = deltaV_Arrival;           % The deltaV @ Arrival [km/s]
                arc.total_delatV     = deltaV_Departure;         % Total deltaV upto this point [km/s] (Summation of all the previous sequences) 
                % JMRM: OLD CODE --> arc.kep        = kep;       % Orbital Keplerian Elements for the trajectory  [a, e, i, Om, om, th] [km, rad]
                arc.kep              = kep_transfer_orbit;       % Orbital Keplerian Elements for the trajectory  [a, e, i, Om, om, th] [km, rad]                                
                arc.departure_r      = departure_r;              % Cartesian Postion @ Departure [km]
                arc.departure_v      = departure_v;              % Cartesian Velocity @ Departure [km/s]
                arc.arrival_r        = arrival_r;                % Cartesian Postion @ Arrival [km]
                arc.arrival_v        = arrival_v;                % Cartesian Velocity @ Arrival [km/s]
                arc.departure_epoch  = curr_departure_epoch;     % Departure Epoch [MJD2000]
                arc.arrival_epoch    = curr_pass_epoch;          % Arrival Epoch [MJD2000]
                arc.tof              = curr_tof;                 % Time of Flight for the Current arc [day]
                arc.total_tof        = curr_tof;                 % Total Time of Flight for the sequence so far [day]
                arc.departure_delbod = curr_departure_celbody;   % Departure Celestial Body for the current arc
                arc.target_celbod    = curr_target_celbody;      % Target Celestial Body for the current arc
                
            else    
                
                %disp(['Creating new ARC '  sprintf('%d', (prev_seq{end}.id + 1))])
                                            
                % Create the new arc 
                arc = TrajectoryArc;
                arc.id               = prev_seq{end}.id + 1;     % Id for the arc as well to control the size of the sequence
                arc.deltaV_departure = deltaV_Departure;         % The deltaV @ Departure [km/s]
                arc.deltaV_arrival   = deltaV_Arrival;           % The deltaV @ Arrival [km/s]
                arc.total_delatV    = prev_seq{end}.total_delatV + deltaV_Departure;% Total deltaV upto this point [km/s] (Summation of all the previous sequences) 
                % JMRM: OLD CODE --> arc.kep        = kep;       % Orbital Keplerian Elements for the trajectory  [a, e, i, Om, om, th] [km, rad]
                arc.kep              = kep_transfer_orbit;       % Orbital Keplerian Elements for the trajectory  [a, e, i, Om, om, th] [km, rad]   
                arc.departure_r      = departure_r;              % Cartesian Postion @ Departure [km]
                arc.departure_v      = departure_v;              % Cartesian Velocity @ Departure [km/s]
                arc.arrival_r        = arrival_r;                % Cartesian Postion @ Arrival [km]
                arc.arrival_v        = arrival_v;                % Cartesian Velocity @ Arrival [km/s]
                arc.departure_epoch  = curr_departure_epoch;     % Departure Epoch [MJD2000]
                arc.arrival_epoch    = curr_pass_epoch;          % Arrival Epoch [MJD2000]
                arc.tof              = curr_tof;                 % Time of Flight for the Current arc [day]
                arc.total_tof        = prev_seq{end}.total_tof + curr_tof;  % Total Time of Flight for the sequence so far [day]
                arc.departure_delbod = curr_departure_celbody;   % Departure Celestial Body for the current arc
                arc.target_celbod    = curr_target_celbody;      % Target Celestial Body for the current arc
            end
            
            % Add the new arc in the current sequence
            curr_seq = prev_seq;
            curr_seq{end + 1} = arc;
             
            % Compute the Minimum ToF required for the next [day]   
            tmp_time = curr_seq{end}.arrival_epoch + global_config.tof_min;
                                     
            % Cotinue the Search Process
            if(    ( ~isempty(curr_missing_celestialbodies) )  ... % Any missing celestial body?    
                && (  tmp_time <= global_config.epoch_end ) )  ... % Enough Time of Flight ? 
                                                    
               % Call recursive again the function
               solutions2 = LambTAN_v1(curr_departure_orbit,         ... % the Next Departure Celestial Body
                                       ind_curr_target_celbody,      ... % Index of the next Departure Celestial Body
                                       curr_seq,                     ... % Current Trajectory Arc
                                       curr_missing_celestialbodies, ... % Missing Celestial Bodies Array                      
                                       global_config,                ... % Global Problem Setting
                                       mu);                          ... % Gravity Parameter              

            end
                
            % Sanity Check
            if length(curr_seq) > (length(global_config.celestialBodies) - 1)
                error('Error found: More number of sequences than elementes to visit!!')
            end
            
            % Add the New sequence into the global Solution Array
            if(length(curr_seq) >= 4)
                solutions{end + 1}  = curr_seq;
            end

            % Increment the Number of found solutions
            num_solutions = num_solutions + 1;
                                    
            % DELETEME
            disp(sprintf('Saving the new Solution sequence #%d with a length of %d', num_solutions, length(curr_seq) + 1 ));                        
                             
        end % End loop for all posible trajectories
        
    end % End loop for Pass Epochs
    
end % End Loop for all missing celestial bodies

end % End of SearchAndCreateNewArc Function

