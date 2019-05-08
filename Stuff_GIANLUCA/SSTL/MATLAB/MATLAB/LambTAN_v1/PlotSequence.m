function [] = PlotSequence(sequence)

% Gravitational parameter [km^2/s^3]
mu = AstroConstants.Sun_Planetary_Const;

% STEP 1: Plot the Earth
% Get Position and Velocity of the DEPARTURE Celestial Body for
% the current departure epoch (firt loop)
curr_celestial = 3;

[departure_r, departure_v] = EphSS( curr_celestial, 0 );

plot_trajectory_v2(departure_r, departure_v, 86400*365, mu, 'g', '', false);

hold on
grid on
%view(120,30)



for i = 1 : length(sequence)

    % Celestial bodies
    celestialBodies = sequence{i}.target_celbod;
                    
    % Departure time
    t_departure = sequence{i}.departure_epoch;

   % Arrival time
    t_arrival = sequence{i}.arrival_epoch;
    
    % Compute the new ToF for the trajectory [days]
    curr_tof = t_arrival - t_departure;
   
%% TRANSFER ORBIT - 1st CASE
%     
%     % This seems to plot the trajectory of the departure celestial body
%     % I thought that solution{1,i}.departure_r and
%     % solution{1,i}.departure_v defined the position and velocity of the
%     % transfer trajectory...
%       plot_trajectory(solution{1,i}.departure_r, solution{1,i}.departure_v,(solution{1,i}.arrival_epoch - solution{1,i}.departure_epoch ) * 86400, mu, 'r'); 
%       
%       
%% TRANSFER ORBIT - 2nd CASE
%     
%     % This seems to plot the trajectory of the departure celestial body
%       plot_trajectory(solution{1,i}.arrival_r, solution{1,i}.arrival_v,(solution{1,i}.arrival_epoch - solution{1,i}.departure_epoch ) * 86400, mu, 'r'); 
%       
%% TRANSFER ORBIT - 3rd CASE
    
    % Get cartesian elements from keplerian elements in solutions{1,i}.kep    
    state = kep2cart(sequence{i}.kep, mu);    
    
    % This plots something that does not goes exactly to the node
    plot_trajectory_v2(state(1:3),state(4:6),-curr_tof * 86400, mu, 'r', sprintf('Dep %d  ', i), true); 
    
%     
%% TRANSFER ORBIT - 4th CASE  
    
%     % Recompute Lambert
% 
%     % Get the Current Departure Celestial Body
%     departure_celbody = celestialBodies{1};
% 
%     % Get the Current Arrival Celestial Body
%     arrival_celbody   = celestialBodies{2};
%  
%     % Set the current departure epoch
%     curr_departure_epoch = t_departure;
% 
%     % Set the current arrival epoch
%     curr_arrival_epoch = t_arrival;
%  
%     % Get Position and Velocity of the DEPARTURE Celestial Body for
%     % the departure epoch 
%     if isa(departure_celbody, 'CelestialBody')
%         [departure_r, departure_v] = StardustTool.CartesianElementsAt( departure_celbody, curr_departure_epoch );
%     elseif isa(departure_celbody, 'double')
%         [departure_r, departure_v] = EphSS( departure_celbody, curr_departure_epoch );
%     end        
% 
%     % Get Position and Velocity of the ARRIVAL Celestial Body for
%     % the arrival epoch 
%     if isa(arrival_celbody, 'CelestialBody')
%         [arrival_r, arrival_v] = StardustTool.CartesianElementsAt( arrival_celbody, curr_arrival_epoch);
%     elseif isa(arrival_celbody, 'double')
%         [arrival_r, arrival_v] = EphSS( arrival_celbody, curr_arrival_epoch );
%     end
% 
% 
%     % Compute the Lambert            
%     [~,~,~,err,vel_initial,vel_final,~,~] = lambertMR(departure_r,      ... % Initial Vector Position
%                                                       arrival_r,        ... % Final position vector
%                                                       curr_tof * 86400, ... % Time of flight [seconds]
%                                                       mu,               ... % Planetary constant of the planet (mu = mass * G) [L^3/T^2]
%                                                       0,                ... % Logical variable defining whether transfer is
%                                                                         ... %   0: direct transfer from R1 to R2 (counterclockwise)
%                                                                         ... %   1: retrograde transfer from R1 to R2 (clockwise)
%                                                       0,                ... % Number of revolutions (Nrev):
%                                                                         ... %   if Nrev = 0 ZERO-REVOLUTION transfer is calculated
%                                                                         ... %   if Nrev > 0 two transfers are possible. Ncase should be
%                                                                         ... %   defined to select one of the two.
%                                                       0,                ... % Logical variable defining the small-a or large-a option in
%                                                                         ... % case of Nrev>0:
%                                                                         ... %   0: small-a option
%                                                                         ... %   1: large-a option
%                                                       0);                   % LambertMR options:
%                                                                             %   optionsLMR(1) = display options:
%                                                                             %     - 0: no display
%                                                                             %     - 1: warnings are displayed only when the algorithm does not converge
%                                                                             %     - 2: full warnings displayed
% 
% 
% 
% 
%  
%         plot_trajectory(departure_r,vel_initial,curr_tof * 86400, mu, 'r'); 


    
%% CELESTIAL BODIES ORBITS
        
    curr_celestial = celestialBodies;

    % Get Position and Velocity of the DEPARTURE Celestial Body for
    % the current departure epoch (firt loop)
    if isa(curr_celestial, 'CelestialBody')
        [departure_r, departure_v] = StardustTool.CartesianElementsAt( curr_celestial, t_departure );
    elseif isa(curr_celestial, 'double')
        [departure_r, departure_v] = EphSS( curr_celestial, t_departure );
    end
    
    %plot_trajectory(departure_r, departure_v, 86400*365, mu, 'b');
    
    plot_trajectory_v2(sequence{i}.arrival_r, sequence{i}.arrival_v, -86400*365, mu, 'b', '', false);

    hold on
    grid on
    %view(120,30)

end


% % STEP 1: Plot the Earth
% % Get Position and Velocity of the DEPARTURE Celestial Body for
% % the current departure epoch (firt loop)
% curr_celestial = 3;
% 
% if isa(curr_celestial, 'CelestialBody')
%     [departure_r, departure_v] = StardustTool.CartesianElementsAt( curr_celestial, t_departure );    
% elseif isa(curr_celestial, 'double')
%     [departure_r, departure_v] = EphSS( curr_celestial, t_departure );
% end
% 
% plot_trajectory_v2(departure_r, departure_v, 86400*365, mu, 'g', '', false);
% 
% hold on
% grid on
% %view(120,30)


