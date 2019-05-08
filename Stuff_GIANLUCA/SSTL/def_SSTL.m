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
function [problem, algo_inner, algo_outer, algo_minmax, algo_decomposition] = def_SSTL()
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
problem.input_minmin_minmax = ''; 
% problem.input_minmin_minmax = 'minmax_IAC1'; 


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
% problem.exact_curves = 3 --> not run ENM
% -------------------------------------------------------------------------
problem.exact_curves = 0;






% -------------------------------------------------------------------------
% load nominal valoue for uncerrtain parameters and optimisa in d
% -------------------------------------------------------------------------
problem.input_uncertain_nominal = [0.0115    0.0536    0.5436    0.0009    2.2160   -1.8571    0.8317    4.4674    0.4200    1.6791    0.1484    0.8327   11.6929    0.5167    4.1687    2.4273  814.4546    0.5481    0.0374    2.7661    2.9267];
 %[0.0100    0.0885    0.6000    0.0010    2.2000    5.0000    0.8000    3.0000    0.5000    1.5000    0.3000    0.8500   10.0000   10.0000   15.0000   50.0000  800.0000    5.0000    5.0000   10.0000   10.0000]; 



% number of sub-functions decomposition
num_functions = 6;  


% number of samples
num_samples = [1 1 1 1 1 5 5 5 5 1 1 1 1 1 1];    % number of samples for each Belief and Plausibility curve of coupled vector



problem.num_functions = num_functions;
for i = 1:problem.num_functions/2*(problem.num_functions-1)
    problem.num_samples{i} = num_samples;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem Definition


%--------------------------------------------------------------------------
% Dimention, Lower and Upper boundaries of the EPISTEMIC parameters 
%--------------------------------------------------------------------------


% dimension (1) satellite, (2) subsystems, (3) mass subsystems, (4) Q_in,
% (5) Q_ex

% delta 
D_axis = 20;
D_e =    0.0012;
D_I =    0.07;
D_RAAN = 30;
D_peri = 0.5;
D_th =   0.025;
D_eff =  0.08;

% avarage values

axis = [68500.3 73250.2 86065.5 49646.4 42049];
e =    [0.902 0.77866 0.513 0.15392 0.001];
I =    [22.81 9.12 1.09 0.36 0.05];
RAAN = [86.63 86.79 85.96 86.85 270];
peri = [180.1 180.6 180.81 180.97 1];
th =   [0 180.08 180.84 4.25 359.95];
eff =  0.8;

problem.lb_u{1} = {[eff-D_eff eff];... % efficiency
              [axis(1)-D_axis axis(1)]; [e(1)-D_e e(1)]; [I(1)-D_I I(1)]; [RAAN(1)-D_RAAN RAAN(1)]; [peri(1)-D_peri peri(1)]; [th(1)-D_th th(1)];... % orbit parameters fire 1
              [axis(2)-D_axis axis(2)]; [e(2)-D_e e(2)]; [I(2)-D_I I(2)]; [RAAN(2)-D_RAAN RAAN(2)]; [peri(2)-D_peri peri(2)]; [th(2)-D_th th(2)];... % orbit parameters fire 2
              [axis(3)-D_axis axis(3)]; [e(3)-D_e e(3)]; [I(3)-D_I I(3)]; [RAAN(3)-D_RAAN RAAN(3)]; [peri(3)-D_peri peri(3)]; [th(3)-D_th th(3)];... % orbit parameters fire 3
              [axis(4)-D_axis axis(4)]; [e(4)-D_e e(4)]; [I(4)-D_I I(4)]; [RAAN(4)-D_RAAN RAAN(4)]; [peri(4)-D_peri peri(4)]; [th(4)-D_th th(4)];... % orbit parameters fire 4
              [axis(5)-D_axis axis(5)]; [e(5)-D_e e(5)]; [I(5)-D_I I(5)]; [RAAN(5)-D_RAAN RAAN(5)]; [peri(5)-D_peri peri(5)]; [th(5)-D_th th(5)];... % orbit parameters fire 5
                         
              }; 

          
problem.ub_u{1} = {[eff eff+D_eff];... % efficiency
              [axis(1) axis(1)+D_axis]; [e(1) e(1)+D_e]; [I(1) I(1)+D_I]; [RAAN(1) RAAN(1)+D_RAAN]; [peri(1) peri(1)+D_peri]; [th(1) th(1)+D_th];... % orbit parameters fire 1
              [axis(2) axis(2)+D_axis]; [e(2) e(2)+D_e]; [I(2) I(2)+D_I]; [RAAN(2) RAAN(2)+D_RAAN]; [peri(2) peri(2)+D_peri]; [th(2) th(2)+D_th];... % orbit parameters fire 2
              [axis(3) axis(3)+D_axis]; [e(3) e(3)+D_e]; [I(3) I(3)+D_I]; [RAAN(3) RAAN(3)+D_RAAN]; [peri(3) peri(3)+D_peri]; [th(3) th(3)+D_th];... % orbit parameters fire 3
              [axis(4) axis(4)+D_axis]; [e(4) e(4)+D_e]; [I(4) I(4)+D_I]; [RAAN(4) RAAN(4)+D_RAAN]; [peri(4) peri(4)+D_peri]; [th(4) th(4)+D_th];... % orbit parameters fire 4
              [axis(5) axis(5)+D_axis]; [e(5) e(5)+D_e]; [I(5) I(5)+D_I]; [RAAN(5) RAAN(5)+D_RAAN]; [peri(5) peri(5)+D_peri]; [th(5) th(5)+D_th];... % orbit parameters fire 5
                           
              };    
          
problem.bpa{1} =  {[.5 .5];...   % efficiency
              [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5];...     
              [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5];...
              [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5];...
              [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5];...
              [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5];...        
              
              };   
          


problem.dim_u = [0 0 0 0 0 1 ...  % u1, u2, u3, u4 ,u5 ,u6
                   0 0 0 0 6 ...  % u12,u13,u14,u15,u16
                     0 0 0 6 ...  % u23,u24,u25,u26            
                       0 0 6 ...  % u34,u35,u36
                         0 6 ...  % u45,u46
                           6 ...  % u56 
            ];    
problem.order_dim_u   = [1];




%--------------------------------------------------------------------------
% Dimention, Lower and Upper boundaries of the DESIGN parameters 
%--------------------------------------------------------------------------

problem.lb_d = [35; 0];     % [D_V, cell] 5000
problem.ub_d = [45.1; 1];   %             10000
 
problem.dim_d = length(problem.lb_d);

%--------------------------------------------------------------------------
% List of the FIXED parameters
%--------------------------------------------------------------------------
problem.fix.day = 59;
problem.fix.MJD0 = 6.939792e+03 -1; %58484+7/24-1;
problem.fix.DoD_max = 0.2;
problem.fix.SoC_max = 0.9;


problem.par_objfun{1}.fix = problem.fix;



%--------------------------------------------------------------------------
% Number of objective functions
%--------------------------------------------------------------------------
problem.n_obj = 1;

%--------------------------------------------------------------------------
% Objective functions;
%--------------------------------------------------------------------------
problem.objfun =     {@sstl_fun};

%--------------------------------------------------------------------------
% Constraints
%--------------------------------------------------------------------------
problem.constraints = {[]}; 

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

