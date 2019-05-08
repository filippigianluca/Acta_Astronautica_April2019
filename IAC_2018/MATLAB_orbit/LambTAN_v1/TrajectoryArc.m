classdef TrajectoryArc
    % TrajectoryArc Summary of this class goes here
    %   Detailed explanation goes here
    
    properties         
        id                  % Id for the arc as well to control the size of the sequence
        deltaV_departure    % The deltaV @ Departure [km/s]
        deltaV_arrival      % The deltaV @ Arrival [km/s]
        total_delatV        % Total deltaV upto this point [km/s] (Summation of all the previous sequences)
        kep                 % Orbital Keplerian Elements for the trajectory  [a, e, i, Om, om, th] [km, rad]
        departure_r         % Cartesian ostion @ Departure [km]
        departure_v         % Cartesian Velocity @ Departure [km/s]
        arrival_r           % Cartesian Postion @ Arrival [km]
        arrival_v           % Cartesian Velocity @ Arrival [km/s]
        departure_epoch     % Departure Epoch [MJD2000]
        arrival_epoch       % Arrival Epoch [MJD2000]
        tof                 % Time of Flight for the Current arc [day]
        total_tof           % Total Time of Flight for the sequence so far [day]
        departure_delbod    % Departure Celestial Body for the current arc
        target_celbod       % Target Celestial Body for the current arc
        lambert_Vel_initial % Initial Lambert Velocity      
        lambert_Vel_final   % Final Lambert Velocity
    end
    
end