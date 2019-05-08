%% GS_meas_2: Ground station measurement model with constant amplitude random noise on spatial components of the spacecraft.
%
% This function provides an estimate of the position of a spacecraft
% according to the coordinates of a given ground station.
% Both input and output states refer to a circular restricted three body
% problem with the Sun and Earth as main bodies.
%
%% Inputs:
%
% * mode : string that defines which part of the function is going to be
%          evaluated. Its values can be:
%          - 'pred' : to be used only in the prediction part of the
%                     navigation filter. Selects the components of the state
%                     that need to be used for the prediction part;
%          - 'event': to be used in the event function only.
%                     According to the time of the integration and the
%                     position of the spacecraft, it only checks whether
%                     the spacecraft is visible from the station or not;
%          - 'meas' : the measurement part. After converting the relative
%                     distance of the spacecraft to the station in the
%                     station's reference frame, it introduces the error
%                     and recomputes the state of the spacecraft in the
%                     original reference frame.
%
% * X_sat_SER : state vector of the spacecraft in the Sun-Earth-Rotating 
%               reference frame.
%
% * t : time of the evaluation, in seconds.
%       It is important in order to evaluate the correct orientation and 
%       position of the station.
%
% * GS : structure containing the parameters of the desired ground station.
%        The correct data structure is provided by ground_station_list.m.
%
% * Earth : structure containing Earth's parameters.
%           The full data structure can be found in the LISA example.
%           The data needed in this function are:
%           - Earth.rot_ang: angles of a 3-1-3 rotation matrix converting
%                            from ECI (Earth Centred Inertial) to SER;
%           - Earth.X      : position vector of the Earth in the SER
%                            reference frame;
%
% * extra : extra contains various parameters useful in LISA navigation and
%           a full description of all the components can be found there.
%           The data needed in this function are:
%           - extra.T_star   : characteristic time of the simulation, in
%                              seconds.
%                              In LISA is the orbital period of the
%                              Sun-Earth-Moon system;
%           - extra.JD_0     : initial Julian Date in seconds;
%           - extra.ell_model: elliptical method used to represent the
%                              Earth. The correct structure is in
%                              ground_station_list.m.
%
% * sensor : structure containing information about the sensors of the
%            spacecraft.
%            The data needed in this function are:
%            - sensor.sensit: sensitivity of the sensor in Km.
%
%% Output:
%
% * out : The output of the function varies according to the mode used.
%         - 'pred'  : position part of the sigma points of the spacecraft;
%         - 'event' : boolean value that is 1 when the spacecraft is
%                     visible from the station and 0 in the other case;
%         - 'meas'  : position vector of the spacecraft.
%
% Author: Francesco Torre
% email: francesco.torre@strath.ac.uk


%% Ground Station measurement model 2: error as constant amplitude noise in X,Y,Z of s/c in SER Rf
function out = GS_meas_2(mode,X_sat_SER,t,GS,Earth,extra,sensor)

switch mode
    
    case 'pred'
        
        H = [1 0 0 0 0 0;
            0 1 0 0 0 0;
            0 0 1 0 0 0];
        
        out = H*X_sat_SER;
        
        
    case 'event'
        
        % Rotation of ECI wrt SER
        rot = -360*t/extra.T_star;

        % SER = Sun-Earth Rotating RF
        R1 = [...
            cosd(Earth.rot_ang(1)+rot)  -sind(Earth.rot_ang(1)+rot)     0;...
            sind(Earth.rot_ang(1)+rot)  cosd(Earth.rot_ang(1)+rot)      0;...
            0                           0                               1];

        R2 = [...
            1   0                       0;...
            0   cosd(Earth.rot_ang(2))  -sind(Earth.rot_ang(2));...
            0   sind(Earth.rot_ang(2))  cosd(Earth.rot_ang(2))];

        R3 = [...
            cosd(Earth.rot_ang(3))  -sind(Earth.rot_ang(3))     0;...
            sind(Earth.rot_ang(3))  cosd(Earth.rot_ang(3))      0;...
            0                       0                           1];

        % Matrix from ECI to SER
        ECI2SER = R3*R2*R1;

        % Position of the satellite in the Earth Centred Inertial (ECI) RF
        X_sat_ECI = ECI2SER.'*(X_sat_SER(1:3,1)-Earth.X);

        % Position of the satellite in the Earth Centred Earth Fixed (ECEF) RF
        R_JD = ECI2ECEF((t+extra.JD_0)/(24*3600));
        X_sat = R_JD*X_sat_ECI;

        % Parameters of Earth model - Using WGS84
        R = extra.ell_model(1).R;
        f = extra.ell_model(1).f;

        % Coordinates of the station in the EBF RF
        % phi = 45;       % deg - latitude
        % lambda = 30;    % deg - longitude
        % h = 10;         % Km  - altitude
        % E0 = 12;        % deg - minimum elevation

        N = R/(sqrt(1-f*(2-f)*sind(GS.phi)^2));

        % Position of the station in the Earth body fixed RF
        X_station = [...
            (N+GS.h)*cosd(GS.phi)*cosd(GS.lambda);...
            (N+GS.h)*cosd(GS.phi)*sind(GS.lambda);...
            ((1-f)^2*N+GS.h)*sind(GS.phi)];

        % Relative position between Earth and the station in the EBF RF
        Rel_pos_EBF = X_sat - X_station;

        % Rotation matrix from EBF to Ground Station (GS)
        E_mat = [...
            -sind(GS.lambda)                  cosd(GS.lambda)                   0;...
            -sind(GS.phi)*cosd(GS.lambda)     -sind(GS.phi)*sind(GS.lambda)     cosd(GS.phi);...
            cosd(GS.phi)*cosd(GS.lambda)      cosd(GS.phi)*sind(GS.lambda)      sind(GS.phi)];

        % Relative position between Earth and the station in the GS RF
        Rel_pos_GS = E_mat*Rel_pos_EBF;

        % Range, Azimuth and Elevation of the satellite in the GS RF
        RAE = [...
            norm(Rel_pos_GS);...
            atan2(Rel_pos_GS(1),Rel_pos_GS(2))*180/pi;...
            atan2(Rel_pos_GS(3),norm(Rel_pos_GS(1:2)))*180/pi];

        % Check if E>E0
        out = RAE(3)-GS.E0 >= 0;
        
        
    case 'meas'
        
        % Rotation of ECI wrt SER
        rot = -360*t/extra.T_star;

        % SER = Sun-Earth Rotating RF
        R1 = [...
            cosd(Earth.rot_ang(1)+rot)  -sind(Earth.rot_ang(1)+rot)     0;...
            sind(Earth.rot_ang(1)+rot)  cosd(Earth.rot_ang(1)+rot)      0;...
            0                           0                               1];

        R2 = [...
            1   0                       0;...
            0   cosd(Earth.rot_ang(2))  -sind(Earth.rot_ang(2));...
            0   sind(Earth.rot_ang(2))  cosd(Earth.rot_ang(2))];

        R3 = [...
            cosd(Earth.rot_ang(3))  -sind(Earth.rot_ang(3))     0;...
            sind(Earth.rot_ang(3))  cosd(Earth.rot_ang(3))      0;...
            0                       0                           1];

        % Matrix from ECI to SER
        ECI2SER = R3*R2*R1;

        % Position of the satellite in the Earth Centred Inertial (ECI) RF
        X_sat_ECI = ECI2SER.'*(X_sat_SER(1:3,1)-Earth.X);

        % Position of the satellite in the Earth Centred Earth Fixed (ECEF) RF
        R_JD = ECI2ECEF((t+extra.JD_0)/(24*3600));
        X_sat = R_JD*X_sat_ECI;

        % Parameters of Earth model
        R = extra.ell_model(1).R;
        f = extra.ell_model(1).f;
        N = R/(sqrt(1-f*(2-f)*sind(GS.phi)^2));

        % Position of the station in the Earth body fixed RF
        X_station = [...
            (N+GS.h)*cosd(GS.phi)*cosd(GS.lambda);...
            (N+GS.h)*cosd(GS.phi)*sind(GS.lambda);...
            ((1-f)^2*N+GS.h)*sind(GS.phi)];

        % Relative position between Earth and the station in the EBF RF
        Rel_pos_EBF = X_sat - X_station;

        % Rotation matrix from EBF to Ground Station (GS)
        E_mat = [...
            -sind(GS.lambda)               cosd(GS.lambda)                0;...
            -sind(GS.phi)*cosd(GS.lambda)     -sind(GS.phi)*sind(GS.lambda)     cosd(GS.phi);...
            cosd(GS.phi)*cosd(GS.lambda)      cosd(GS.phi)*sind(GS.lambda)      sind(GS.phi)];

        % Relative position between Earth and the station in the GS RF
        Rel_pos_GS = E_mat*Rel_pos_EBF;

        % Range, Azimuth and Elevation of the satellite in the GS RF
        RAE = [...
            norm(Rel_pos_GS);...
            atan2(Rel_pos_GS(1),Rel_pos_GS(2))*180/pi;...
            atan2(Rel_pos_GS(3),norm(Rel_pos_GS(1:2)))*180/pi];
        
        if RAE(3) < GS.E0
            out = NaN(3,1);
            disp('Something bad happened')
            return;
        end
        
        % Noise: random error in the three components with fixed amplitude
        noise = sensor.sensit*rand(3,1).*sign(randn(3,1));
        
        % Converting the noised measure in the GS RF
        Rel_pos_GS_n = Rel_pos_GS + noise;

        % Converting the noised measure in the EBF RF
        Rel_pos_EBF_n = (E_mat.')*Rel_pos_GS_n;

        % Noised position of the satellite in the ECEF RF
        X_sat_n = Rel_pos_EBF_n + X_station;

        % Converting the noised measure in the ECI RF
        X_sat_ECI_n = (R_JD.')*X_sat_n;

        % Converting the noised measure in the SER RF
        out = ECI2SER*X_sat_ECI_n + Earth.X;

end