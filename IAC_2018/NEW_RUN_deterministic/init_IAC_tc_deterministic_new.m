% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of run of a min-max problem using MP_AIDEA in a single objective
% optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [problem] = init_IAC_tc_deterministic_new()
% "def_ebro_problem_cubesat" for PrmaDiNatale folder



%% Type of Input
problem.flag_input.dmin_load = 0;             % load (1) or not (0) the design vector for best-case problem         
problem.flag_input.dmax_load = 1;             % load (1) or not (0) the design vector for worst-case problem        
problem.flag_input.umin_load = 0;             % load (1) or not (0) the uncertain vector for best-case problem   
problem.flag_input.umax_load = 1;             % load (1) or not (0) the uncertain vector for worst-case problem   
problem.flag_input.fmin_load = 0;             % load (1) or not (0) the value of objective function f(dmin_load, umin_load) 
problem.flag_input.fmax_load = 0;             % load (1) or not (0) the value of objective function f(dmax_load, umax_load) 




%% Type of Output
problem.flag_output.u_run    = 0;              % Evaluate (1) or not (0) the maximum of f for a fixed design d
problem.flag_output.dmin_run = 0;              % Evaluate (1) or not (0) the minimum of f for a fixed design u
problem.flag_output.minmin   = 0;              % Evaluate (1) or not (0) the Best Case Scenario (min-min problem)
problem.flag_output.minmax   = 0;              % Evaluate (1) or not (0) the Worst Case Scenario (min-max problem)
problem.flag_output.Belief   = 1;              % Use (1) or not (0) Evidence-Network-Model to reconstruct the Belief curve
problem.flag_output.Plausibility = 0;          % Use (1) or not (0) Evidence-Network-Model to reconstruct the Plausibility curve
problem.flag_output.exact_Belief = 0;          % Evaluate (1) or not (0) the exact Belief curve for the fixed design d
problem.flag_output.exact_Plausibility = 0;    % Evaluate (1) or not (0) the exact Plausibility curve for the fixed design d
problem.flag_output.plot = 1;                  % plot (1) or not (0) the curves




%% Load the inputs
problem.input.d_min    = [];    
problem.input.u_min{1} = [];  
% MAX MASS
% load('d_minmax_nu400.mat','dminmax')
% load('u_minmax_nu400.mat','uminmax')
load('minmax_d_nu600_13_04','minmax_d');
load('minmax_u_nu600_13_04','minmax_u');
problem.input.d_max    = minmax_d;
problem.input.u_max{1} = minmax_u;
problem.input.F_max{1} = [];

% % MIN DATA VOLUME
% problem.input.u_max{1} = [
% problem.input.F_max{1} = []; 

load('u_max_C_nu400.mat','umaxC')
problem.input.u_sample{1} = []; 

problem.input.F_min{1} = [];    
load('u_nominal.mat','u_nominal')
problem.input.vector_u_nominal = [u_nominal(1:12) u_nominal(14:end)];





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem Definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% objectives and constraints
problem.n_obj      = 1;
problem.objfun     = {@IAC2018_only_mass2_for_belief};            %each function is in the form [f_i] = objfun(i)(d,u,par)
problem.constraint = {[]}; 
problem.obj_constr = 0; 


%% design variables
problem.lb_d = [10 30  ...     % AOCS 
                 7  0 0 ...    % TTC  
                 0  3 0 1 ...  % EPS
                 1 0 ...    % PAYLOAD
                 0 ]';         % OBDH

problem.ub_d = [60 90 ...      % AOCS
                10  1 1 ...    % TTC 
                 1  5 1 3 ...  % EPS
                 5 1 ...
                 1 ]';
 
problem.dim_d = length(problem.lb_d);



%% uncertain variables
problem.dim_u   = 20;
problem.lb_u{1} = { ... 
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
%                     [0 10]; ...               %       * ep(7)  = uncertainty on array temperature, C                    (ALICINO uncertain, [0 10][10 15])
                    [0 10]; ...               %       * ep(10) = uncertainty on power requirements, %                   (ALICINO uncertain, [0 10][10 20])
                    [0 15]; ...               %       * ep(13) = uncertainty on Delta rho_sa [% kg/m^2] density solar array, [0 15][15 30]
                    [0 50]; ...               %       * ep(14) = uncertainty on Delta D_cell [%] degradation solar array,    [0 50][50 100]

                    [600 800]; ...            % PAYLOAD * altitude
                    [0 5]; ...                %         * elevation angle
                    [0 5]; ...                %         * inclination (%)
                    
                    [0 10]; ...               % OBDH  * mass (%)
                    [0 10]; ...               %       * power (%)
                  };
 
 
                                              
problem.ub_u{1} = { ... 
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
%                     [10 15]; ...              %       * ep(7)  = uncertainty on array temperature, C                    (ALICINO uncertain, [0 10][10 15])
                    [10 20]; ...              %       * ep(10) = uncertainty on power requirements, %                   (ALICINO uncertain, [0 10][10 20])
                    [15 30]; ...              %       * ep(13) = uncertainty on Delta rho_sa [% kg/m^2] density solar array, [0 15][15 30]
                    [50 100]; ...             %       * ep(14) = uncertainty on Delta D_cell [%] degradation solar array,    [0 50][50 100]    

                    [800 1000]; ...           % PAYLOAD
                    [5 10]; ...
                    [5 10]; ...
                    
                    [10 20]; ...              % OBDH
                    [10 20]; ...
                  };     

problem.bpa{1}  = {[0.5 0.5]; ...             % AOCS
                   [0.5 0.5]; ...
                   [0.5 0.5]; ...
                   [0.5 0.5]; ...
                   [0.4 0.6]; ...
                   [0.5 0.5]; ...
                   
                   [0.3 0.7]; ...            % TTC
                   [0.3 0.7]; ...
                   [0.3 0.7]; ...
                   [0.4 0.6]; ...
                   [0.4 0.6]; ...
                   
                   [0.4 0.6]; ...            % POWER
%                    [0.4 0.6]; ...
                   [0.3 0.7]; ...
                   [0.5 0.5]; ...
                   [0.8 0.2]; ...
                   
                   [0.4 0.6]; ...            % PAYLOAD
                   [0.4 0.6]; ...
                   [0.3 0.7]; ...
                   
                   [0.5 0.5]; ...            % OBDH
                   [0.8 0.2]; ...
                  };


%% fixed variables
% problem.par_objfun{1}.fix.prova = 2;
par_objfun.fix.time = 365;
par_objfun.fix.nu = 600;
    problem.par_objfun = {par_objfun};



%% network data
problem.num_functions = 5; % number of sub-functions decomposition

% coupling matrix
%                      f1 f2
%
% problem.CM =  u1    [1  0; ...
%               u2     0  1; ...
%               u3     1  1]

problem.CM = [1 0 1 0 0; ... % 1
              1 0 1 0 0; ... % 2
              1 0 0 0 0; ... % 3
              1 0 0 0 0; ... % 4
              1 0 1 0 0; ... % 5
              1 0 1 0 0; ... % 6
              0 1 0 0 0; ... % 7
              0 1 1 0 0; ... % 8
              0 1 1 0 0; ... % 9
              0 1 1 0 0; ... % 10
              0 1 0 0 0; ... % 11
              0 0 1 0 0; ... % 12
              0 0 1 0 0; ... % 14
              0 0 1 0 0; ... % 15
              0 0 1 0 0; ... % 16
              1 0 0 1 0; ... % 17
              1 0 0 1 0; ... % 18
              1 0 0 1 0; ... % 19
              0 0 0 0 1; ... % 20
              0 0 1 0 1];    % 21
             
          
problem.dim_u_i       = [2 2 4 0 1 0 4 3 0 3 0 0 0 1];     % [u1, u2, u3, u12, u13, u23]
problem.order_dim_u   = [3 4 ...              % 1         AOCS
                         7 11 ...             % 2         TTC
                         12 13 14 15 ...   % 3         POWER
                         ...                  % 4         PAYLOAD
                         19 ...               % 5         OBDH
                         ...                  % 1-2
                         1 2 5 6 ...          % 1-3
                         16 17 18 ...         % 1-4
                         ...                  % 1-5
                         8 9 10 ...           % 2-3
                         ...                  % 2-4
                         ...                  % 2-5
                         ...                  % 3-4
                         20 ...               % 3-5
                         ];


% number of samples
num_samples = 3;
for i = 1:problem.num_functions/2*(problem.num_functions-1)
    problem.num_samples{i} = num_samples;
end





return

