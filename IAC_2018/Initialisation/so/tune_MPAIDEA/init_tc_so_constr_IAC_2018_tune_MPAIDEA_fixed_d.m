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
function problem = init_tc_so_constr_IAC_2018_tune_MPAIDEA_fixed_d()
% "def_ebro_problem_cubesat" for PrmaDiNatale folder



global nfevalglobal;
nfevalglobal = 0;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose the output

% -------------------------------------------------------------------------
% problem.output = 0 --> only minmax-Belief
% problem.output = 1 --> only minmin-Plausibility
% problem.output = 2 --> both minmax-Belief and minmin-Plausibility
% -------------------------------------------------------------------------
problem.output = 0;


% -------------------------------------------------------------------------
% problem.input = 0 --> do minmax and minmin
% problem.input = 1 --> load d, u_min, u_max
% problem.input = 2 --> load d, run max and min
% problem.input = 3 --> load u, run min d
% -------------------------------------------------------------------------
problem.input = 0;


% -------------------------------------------------------------------------
% if input == 1 or input == 2.
% if you have the structure(s) 'minmax' or 'minmin' from the 
% optimise_macs_minmax problem write the name of the file .mat here:
% -------------------------------------------------------------------------
% problem.input_minmin_minmax = ''; 
problem.input_minmin_minmax = ''; 


% -------------------------------------------------------------------------
% oterwise problem.input.problem.input.minmin_minmax = [] and:
% -------------------------------------------------------------------------
problem.input_d_min    = [];    % d_min    = [d1, d2, ..., dn] 
problem.input_d_max    = [];    % d_max    = [d1, d2, ..., dn] (lower bound)
problem.input_u_min    = [];    % u_min{1} = [u1, u2, ...,un]
problem.input_u_max{1} = [];
problem.input_F_min    = [];    % F_min    = F 
problem.input_F_max    = [];    % F_max    = F



% -------------------------------------------------------------------------
% do decompoition AND/OR exact Belief OR run max_u
% problem.exact_curves = 0 --> NO exact curves
% problem.exact_curves = 1 --> YES exact curves
% problem.exact_curves = 2 --> ONLY exact curves
% problem.exact_curves = 3 --> DO NOT run decompoition (ONLY max_u for a
% given d)
% -------------------------------------------------------------------------
problem.exact_curves = 3;





% -------------------------------------------------------------------------
% load nominal valoue for uncerrtain parameters and optimisa in d
% -------------------------------------------------------------------------
problem.input_uncertain_nominal = [];




% number of sub-functions decomposition
num_functions = 5;  


% number of samples
num_samples = 1;    




problem.num_functions = num_functions;
for i = 1:problem.num_functions/2*(problem.num_functions-1)
    problem.num_samples{i} = num_samples;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem Definition


%--------------------------------------------------------------------------
% Dimention, Lower and Upper boundaries of the EPISTEMIC parameters 
%--------------------------------------------------------------------------


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
                    [0 10]; ...               %       * ep(7)  = uncertainty on array temperature, C                    (ALICINO uncertain, [0 10][10 15])
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
                    [10 15]; ...              %       * ep(7)  = uncertainty on array temperature, C                    (ALICINO uncertain, [0 10][10 15])
                    [10 20]; ...              %       * ep(10) = uncertainty on power requirements, %                   (ALICINO uncertain, [0 10][10 20])
                    [15 30]; ...              %       * ep(13) = uncertainty on Delta rho_sa [% kg/m^2] density solar array, [0 15][15 30]
                    [50 100]; ...             %       * ep(14) = uncertainty on Delta D_cell [%] degradation solar array,    [0 50][50 100]    

                    [800 1000]; ...           % PAYLOAD
                    [5 10]; ...
                    [5 10]; ...
                    
                    [10 20]; ...              % OBDH
                    [10 20]; ...
                  };     

problem.bpa{1}  = {[0.5 0.5]; ...
                   [0.5 0.5]; ...
                   [0.5 0.5]; ...
                   [0.5 0.5]; ...
                   [0.4 0.6]; ...
                   [0.5 0.5]; ...
                   [0.3 0.7]; ...
                   [0.3 0.7]; ...
                   [0.3 0.7]; ...
                   [0.4 0.6]; ...
                   [0.4 0.6]; ...
                   [0.4 0.6]; ...
                   [0.4 0.6]; ...
                   [0.3 0.7]; ...
                   [0.5 0.5]; ...
                   [0.8 0.2]; ...
                   [0.4 0.6]; ...
                   [0.4 0.6]; ...
                   [0.3 0.7]; ...
                   [0.5 0.5]; ...
                   [0.8 0.2]; ...
                  };

problem.dim_u         = [2 2 5 0 1 0 4 3 0 3 0 0 0 1];     % [u1, u2, u3, u12, u13, u23]
problem.order_dim_u   = [3 4 ...              % 1
                         7 11 ...             % 2
                         12 13 14 15 16 ...   % 3
                         ...                  % 4
                         20 ...            % 5
                         ...                  % 1-2
                         1 2 5 6 ...          % 1-3
                         17 18 19 ...         % 1-4
                         ...                  % 1-5
                         8 9 10 ...           % 2-3
                         ...                  % 2-4
                         ...                  % 2-5
                         ...                  % 3-4
                         21 ...            % 3-5
                         ];




%--------------------------------------------------------------------------
% Dimention, Lower and Upper boundaries of the DESIGN parameters 
%--------------------------------------------------------------------------

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

%--------------------------------------------------------------------------
% List of the FIXED parameters
%--------------------------------------------------------------------------
problem.fix = [];
problem.fix.time = 365;
problem.fix.nu = 400;


problem.par_objfun.fix = problem.fix;



%--------------------------------------------------------------------------
% Number of objective functions
%--------------------------------------------------------------------------
problem.n_obj = 1;

%--------------------------------------------------------------------------
% Objective functions;
%--------------------------------------------------------------------------
problem.objfun =     {@IAC_function_MASS_tune_fix_d};

%--------------------------------------------------------------------------
% Constraints
%--------------------------------------------------------------------------
problem.constraints = {@IAC_constr_tune_fix_d}; 

problem.func_constraints = 0;






% design variables
problem.dim_d = problem.dim_d;
problem.lb_d = problem.lb_d;
problem.ub_d = problem.ub_d;
% uncertain variables
problem.dim_u_i = problem.dim_u;
dim_u = sum(problem.dim_u_i);
problem.dim_u = dim_u;

for n =1:problem.n_obj
    problem.lb_u{n} = problem.lb_u{n};
    problem.ub_u{n} = problem.ub_u{n};
end









return

