function [lb_d, ub_d, lb_u, ub_u] = bound_AOCS()

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

lb_d = [ 15, ...  %       * x(1) = pointing accuracy, deg 
         0, ...   %       * x(2) = offset between centre of gravity and centre of pressure, m
         0.034, ...   %   * x(3) = area normal to velocity vector, m^2  (ALICINO)
         0.034, ...   %   * x(4) = area normal to sun line, m^2         (ALICINO)
         0.5, ... %       * x(5) = reflectance factor                   (ALICINO)
         0.0005,  ...  %  * x(6) = spacecraft residual dipole, A*m^2    (ALICINO)
         2, ...   %       * x(7) = drag coefficient                     (ALICINO)
         0.7, ...   %     * x(8) = actuator type: [0 0.5) = magtorq, [0.5 1] = thruster
         100, ...   %     * x(9) = specific impulse, s
         10, ...   %      * x(10) = burn time for momentum dumping, s
         10, ...  %       * x(11) = slew angle, deg
         30, ...  %       * x(12) = time for slew maneuvre, s
         10, ...   %      * x(13) = initial tumbling rate, rad/s
         -100 ...    %    * x(14) = spin angle, deg 
];

     
     
ub_d = [40, ...   %        * x(1) = pointing accuracy, deg
         0.3, ... %        * x(2) = offset between centre of gravity and centre of pressure, m 
         0.15, ...   %        * x(3) = area normal to velocity vector, m^2
         0.15, ...   %        * x(4) = area normal to sun line, m^2
         0.7, ...   %      * x(5) = reflectance factor
         0.0015, ...   %   * x(6) = spacecraft residual dipole, A*m^2
         2.5, ...  %       * x(7) = drag coefficient 
         0.8, ...  %       * x(8) = actuator type: [0 0.5) = magtorq, [0.5 1] = thruster
         200, ...  %       * x(9) = specific impulse, s
         50, ...  %       * x(10) = burn time for momentum dumping, s
         60, ...   %       * x(11) = slew angle, deg
         90, ...   %       * x(12) = time for slew maneuvre, s
         20, ...   %       * x(13) = initial tumbling rate, rad/s
         100 ...    %       * x(14) = spin angle, deg 
];

lb_u = [5500, ...       %       * ep(1) = orbit period, s
        3e-6, ...       %       * ep(2) = mean planet's magnetic field stregth, T
        500, ...        %       * ep(3) = mean incident solar radiation, W/m^2
        1e-9, ...       %       * ep(4) = mean dynamic pressure, Pa
        1e-7, ...       %       * ep(5) = mean gravitational field strength, deg
        0, ...          %       * ep(6) = angle between spacecraft z axis and nadir vector, deg
        0.04, ...          %       * ep(7) = spacecraft moment of inertia about x axis, kg*m^2
        0.01, ...          %       * ep(8) = spacecraft moment of inertia about y axis, kg*m^2
        0.01, ...          %       * ep(9) = spacecraft moment of inertia about z axis, kg*m^2
        1, ...         %       * ep(10) = thruster moment arm, m
        1, ...         %       * ep(11) = number of thrusters
        1, ...         %       * ep(12) = number of dumping maneuvres during S/C lifetime
        1, ...         %       * ep(13) = number of slew maneuvres during S/C lifetime
];



ub_u = [6500, ...      %       * ep(1) = orbit period, s
        3e-3, ...      %       * ep(2) = mean planet's magnetic field stregth, T
        2500, ...      %       * ep(3) = mean incident solar radiation, W/m^2
        1e-7, ...        %       * ep(4) = mean dynamic pressure, Pa
        1e-4, ...       %       * ep(5) = mean gravitational field strength, deg
        20, ...       %       * ep(6) = angle between spacecraft z axis and nadir vector, deg
        1, ...         %       * ep(7) = spacecraft moment of inertia about x axis, kg*m^2
        1, ...         %       * ep(8) = spacecraft moment of inertia about y axis, kg*m^2
        1, ...         %       * ep(9) = spacecraft moment of inertia about z axis, kg*m^2
        2, ...         %       * ep(10) = thruster moment arm, m
        5, ...         %       * ep(11) = number of thrusters
        5, ...         %       * ep(12) = number of dumping maneuvres during S/C lifetime
        5, ...         %       * ep(13) = number of slew maneuvres during S/C lifetime
];

return