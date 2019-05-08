% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [problem_minmax] = init_tc_mo_IAC2018(varargin)

if (nargin == 0)
    id = 1;
else
    id = varargin{1};
end



%% initialise problem TC1

%--------------------------------------------------------------------------
% type of problem
%--------------------------------------------------------------------------
problem_minmax.sign_inner = 1;          % -1 will run minmin
problem_minmax.flag_output.minmin = 0;    
problem_minmax.flag_output.minmax = 1;    
problem_minmax.obj_constr = 0;  
%--------------------------------------------------------------------------
% objectives
%--------------------------------------------------------------------------
problem_minmax.n_obj = 2;
problem_minmax.objfun = {@CUBESAT_5subsystems_MASS, @CUBESAT_5subsystems_RES};            %each function is in the form [f_i] = objfun(i)(d,u,par)
problem_minmax.constraints = {[], []};
problem_minmax.par_objfun = {struct, struct};


%--------------------------------------------------------------------------
% design variables
%--------------------------------------------------------------------------
dim_d = 12;
problem_minmax.dim_d = dim_d;
problem_minmax.lb_d = [10 30  ...     % AOCS 
                         7  0 0 ...    % TTC  
                         0  3 0 1 ...  % EPS
                         1 0 ...    % PAYLOAD
                         0 ]';         % OBDH
problem_minmax.ub_d =  [60 90 ...      % AOCS
                        10  1 1 ...    % TTC 
                         1  5 1 3 ...  % EPS
                         5 1 ...
                         1 ]';

                     
%--------------------------------------------------------------------------
% uncertain variables
%--------------------------------------------------------------------------
dim_u = 21;
problem_minmax.dim_u = dim_u;
problem_minmax.lb_u = {{ ... 
                    [0.005 0.01]; ...         % AOCS  * x(2) = offset between centre of gravity and centre of pressure, m   (ALICINO uncertain, [0.005 0.01][0.01 0.02])
                    [0.034 0.0885]; ...       %       * x(3) = area normal to velocity vector, m^2                          (ALICINO uncertain, [0.034 0.0885][0.0885 0.15])
                    [0.5 0.6]; ...            %       * x(5) = reflectance factor                                           (ALICINO uncertain, [0.5 0.6][0.6 0.7])
                    [0.0005 0.001]; ...       %       * x(6) = spacecraft residual dipole, A*m^2                            (ALICINO uncertain, [0.0005 0.001][0.001 0.0015])
                    [2 2.2]; ...              %       * x(7) = drag coefficient                                             (ALICINO uncertain, [2 2.2][2.2 2.5])                       
                    [-10 5]; ...              %       * ep(14) = Delta Inertia                                              (ALICINO uncertain, [-10 5][5 10])                       
                    
                    [0.6 0.8]; ...            % TTC   * x(3) = antenna efficiency                       (ALICINO uncertainty  [0.6 0.8][0.8 0.9])
                    [1 3]; ...                %       * x(4) = antenna gain, dB                         (ALICINO uncertainty  [1 3][3 5] trasmit AG)
                    [0.1 0.5]; ...            %       * x(5) = onboard loss, dB                         (ALICINO uncertainty  [0.1 0.5][0.5 1])
                    [0.5 1.5]; ...            %       * x(6) = other unmodelled losses, dB              (ALICINO uncertainty  [0.5 1.5][1.5 2])
                    [0.1 0.3]; ...            %       * x(7) = mass of distribution network, kg         (ALICINO uncertainty  [0.1 0.3][0.2 0.5])
                    
                    [0.8 0.85]; ...           % POWER * x(4)   = cell packing efficiency (or assembly factor) [0,1]     (ALICINO uncertain, [0.8 0.85][0.85 0.9])   
                    [0 10]; ...               %       * ep(7)  = uncertainty on array temperature, C                    (ALICINO uncertain, [0 10][10 15])
                    [0 10]; ...               %       * ep(10) = uncertainty on power requirements, %                   (ALICINO uncertain, [0 10][10 20])
                    [0 15]; ...               %       * ep(13) = uncertainty on Delta rho_sa [% kg/m^2] density solar array, [0 15][15 30]
                    [0 50]; ...               %       * ep(14) = uncertainty on Delta D_cell [%] degradation solar array,    [0 50][50 100]

                    [600 800]; ...            % PAYLOAD * altitude
                    [0 5]; ...                %         * elevation angle
                    [0 5]; ...                %         * inclination (%)
                    
                    [0 10]; ...               % OBDH  * mass (%)
                    [0 10]; ...               %       * power (%)
                  }, ...
                  { ... 
                    [0.005 0.01]; ...         % AOCS  * x(2) = offset between centre of gravity and centre of pressure, m   (ALICINO uncertain, [0.005 0.01][0.01 0.02])
                    [0.034 0.0885]; ...       %       * x(3) = area normal to velocity vector, m^2                          (ALICINO uncertain, [0.034 0.0885][0.0885 0.15])
                    [0.5 0.6]; ...            %       * x(5) = reflectance factor                                           (ALICINO uncertain, [0.5 0.6][0.6 0.7])
                    [0.0005 0.001]; ...       %       * x(6) = spacecraft residual dipole, A*m^2                            (ALICINO uncertain, [0.0005 0.001][0.001 0.0015])
                    [2 2.2]; ...              %       * x(7) = drag coefficient                                             (ALICINO uncertain, [2 2.2][2.2 2.5])                       
                    [-10 5]; ...              %       * ep(14) = Delta Inertia                                              (ALICINO uncertain, [-10 5][5 10])                       
                    
                    [0.6 0.8]; ...            % TTC   * x(3) = antenna efficiency                       (ALICINO uncertainty  [0.6 0.8][0.8 0.9])
                    [1 3]; ...                %       * x(4) = antenna gain, dB                         (ALICINO uncertainty  [1 3][3 5] trasmit AG)
                    [0.1 0.5]; ...            %       * x(5) = onboard loss, dB                         (ALICINO uncertainty  [0.1 0.5][0.5 1])
                    [0.5 1.5]; ...            %       * x(6) = other unmodelled losses, dB              (ALICINO uncertainty  [0.5 1.5][1.5 2])
                    [0.1 0.3]; ...            %       * x(7) = mass of distribution network, kg         (ALICINO uncertainty  [0.1 0.3][0.2 0.5])
                    
                    [0.8 0.85]; ...           % POWER * x(4)   = cell packing efficiency (or assembly factor) [0,1]     (ALICINO uncertain, [0.8 0.85][0.85 0.9])   
                    [0 10]; ...               %       * ep(7)  = uncertainty on array temperature, C                    (ALICINO uncertain, [0 10][10 15])
                    [0 10]; ...               %       * ep(10) = uncertainty on power requirements, %                   (ALICINO uncertain, [0 10][10 20])
                    [0 15]; ...               %       * ep(13) = uncertainty on Delta rho_sa [% kg/m^2] density solar array, [0 15][15 30]
                    [0 50]; ...               %       * ep(14) = uncertainty on Delta D_cell [%] degradation solar array,    [0 50][50 100]

                    [600 800]; ...            % PAYLOAD * altitude
                    [0 5]; ...                %         * elevation angle
                    [0 5]; ...                %         * inclination (%)
                    
                    [0 10]; ...               % OBDH  * mass (%)
                    [0 10]; ...               %       * power (%)
                  }};
problem_minmax.ub_u = {{ ... 
                    [0.01 0.02]; ...          % AOCS  * x(2) = offset between centre of gravity and centre of pressure, m   (ALICINO uncertain, [0.005 0.01][0.01 0.02])
                    [0.0885 0.15]; ...        %       * x(3) = area normal to velocity vector, m^2                          (ALICINO uncertain, [0.034 0.0885][0.0885 0.15])
                    [0.6 0.7]; ...            %       * x(5) = reflectance factor                                           (ALICINO uncertain, [0.5 0.6][0.6 0.7])
                    [0.001 0.0015]; ...       %       * x(6) = spacecraft residual dipole, A*m^2                            (ALICINO uncertain, [0.0005 0.001][0.001 0.0015])
                    [2.2 2.5]; ...            %       * x(7) = drag coefficient                                             (ALICINO uncertain, [2 2.2][2.2 2.5])                       
                    [5 10];...                %       * ep(14) = Delta Inertia                                              (ALICINO uncertain, [-10 5][5 10])     
                    
                    [0.8 0.9];...             % TTC   * x(3) = antenna efficiency                       (ALICINO uncertainty  [0.6 0.8][0.8 0.9])
                    [3 5]; ...                %       * x(4) = antenna gain, dB                         (ALICINO uncertainty  [1 3][3 5] trasmit AG)
                    [0.5 1]; ...              %       * x(5) = onboard loss, dB                         (ALICINO uncertainty  [0.1 0.5][0.5 1])
                    [1.5 2]; ...              %       * x(6) = other unmodelled losses, dB              (ALICINO uncertainty  [0.5 1.5][1.5 2])
                    [0.2 0.5]; ...            %       * x(7) = mass of distribution network, kg         (ALICINO uncertainty  [0.1 0.3][0.2 0.5])
                    
                    [0.85 0.9]; ...           % POWER * x(4)   = cell packing efficiency (or assembly factor) [0,1]     (ALICINO uncertain, [0.8 0.85][0.85 0.9])   
                    [10 15]; ...              %       * ep(7)  = uncertainty on array temperature, C                    (ALICINO uncertain, [0 10][10 15])
                    [10 20]; ...              %       * ep(10) = uncertainty on power requirements, %                   (ALICINO uncertain, [0 10][10 20])
                    [15 30]; ...              %       * ep(13) = uncertainty on Delta rho_sa [% kg/m^2] density solar array, [0 15][15 30]
                    [50 100]; ...             %       * ep(14) = uncertainty on Delta D_cell [%] degradation solar array,    [0 50][50 100]    

                    [800 1000]; ...           % PAYLOAD
                    [5 10]; ...
                    [5 10]; ...
                    
                    [10 20]; ...              % OBDH
                    [10 20]; ...
                  }, ...
                  { ... 
                    [0.01 0.02]; ...          % AOCS  * x(2) = offset between centre of gravity and centre of pressure, m   (ALICINO uncertain, [0.005 0.01][0.01 0.02])
                    [0.0885 0.15]; ...        %       * x(3) = area normal to velocity vector, m^2                          (ALICINO uncertain, [0.034 0.0885][0.0885 0.15])
                    [0.6 0.7]; ...            %       * x(5) = reflectance factor                                           (ALICINO uncertain, [0.5 0.6][0.6 0.7])
                    [0.001 0.0015]; ...       %       * x(6) = spacecraft residual dipole, A*m^2                            (ALICINO uncertain, [0.0005 0.001][0.001 0.0015])
                    [2.2 2.5]; ...            %       * x(7) = drag coefficient                                             (ALICINO uncertain, [2 2.2][2.2 2.5])                       
                    [5 10];...                %       * ep(14) = Delta Inertia                                              (ALICINO uncertain, [-10 5][5 10])     
                    
                    [0.8 0.9];...             % TTC   * x(3) = antenna efficiency                       (ALICINO uncertainty  [0.6 0.8][0.8 0.9])
                    [3 5]; ...                %       * x(4) = antenna gain, dB                         (ALICINO uncertainty  [1 3][3 5] trasmit AG)
                    [0.5 1]; ...              %       * x(5) = onboard loss, dB                         (ALICINO uncertainty  [0.1 0.5][0.5 1])
                    [1.5 2]; ...              %       * x(6) = other unmodelled losses, dB              (ALICINO uncertainty  [0.5 1.5][1.5 2])
                    [0.2 0.5]; ...            %       * x(7) = mass of distribution network, kg         (ALICINO uncertainty  [0.1 0.3][0.2 0.5])
                    
                    [0.85 0.9]; ...           % POWER * x(4)   = cell packing efficiency (or assembly factor) [0,1]     (ALICINO uncertain, [0.8 0.85][0.85 0.9])   
                    [10 15]; ...              %       * ep(7)  = uncertainty on array temperature, C                    (ALICINO uncertain, [0 10][10 15])
                    [10 20]; ...              %       * ep(10) = uncertainty on power requirements, %                   (ALICINO uncertain, [0 10][10 20])
                    [15 30]; ...              %       * ep(13) = uncertainty on Delta rho_sa [% kg/m^2] density solar array, [0 15][15 30]
                    [50 100]; ...             %       * ep(14) = uncertainty on Delta D_cell [%] degradation solar array,    [0 50][50 100]    

                    [800 1000]; ...           % PAYLOAD
                    [5 10]; ...
                    [5 10]; ...
                    
                    [10 20]; ...              % OBDH
                    [10 20]; ...
                  }};

%--------------------------------------------------------------------------
% maxnfeval: optional to give it here or in the algorithm
%--------------------------------------------------------------------------
% problem_minmax.maxnfeval = 2.0e1;



%--------------------------------------------------------------------------
% List of the FIXED parameters
%--------------------------------------------------------------------------
problem_minmax.fix = [];
problem_minmax.fix.time = 365;


problem_minmax.par_objfun{1}.fix = problem_minmax.fix;
problem_minmax.par_objfun{2}.fix = problem_minmax.fix;

return