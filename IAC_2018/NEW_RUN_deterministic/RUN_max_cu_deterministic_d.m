% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ----------
%--------------- e-mail: smart@strath.ac.uk -------------------------------
%------------------- Authors: SMART developers team -----------------------
% clear all; close all; clc



warning('off','all')
%warning('on','all')

%% ------------------------------------------------------------------------
% add folders to path
%--------------------------------------------------------------------------

repo_folder = '../../';

addpath(genpath([repo_folder 'smart-o2c/Optimisation']))
addpath(genpath([repo_folder 'smart-o2c/Other_Tools']))
addpath(genpath([repo_folder 'smart-o2c/Problems/EBRO']))
rmpath(genpath([repo_folder 'smart-o2c/Problems/EBRO/CEC2019']))
addpath(genpath([repo_folder 'utopiae-reliables/IAC_2018']))
addpath(genpath([repo_folder 'utopiae-reliables/Stuff_DANDA']))
rmpath(genpath([repo_folder 'utopiae-reliables/IAC_2018/Optimisation_IAC2018']))

%--------------------------------------------------------------------------



% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

global nfevalglobal;
nfevalglobal = 0;




%% initialisation
init_problem   = str2func('init_IAC_tc_ebro');
init_algorithm = str2func('init_IAC_algo_so_ebro');

savefolder = strcat('RESULTS/SO/');

save_flag = 1;


disp(strcat('ebro_run_'));


%% initialise problem
[ problem_ebro ] = init_problem();
%% ------------------------------------------------------------------------
% Type of Input
problem_ebro.flag_input.dmin_load = 0;             % load (1) or not (0) the design vector for best-case problem         
problem_ebro.flag_input.dmax_load = 1;             % load (1) or not (0) the design vector for worst-case problem        
problem_ebro.flag_input.umin_load = 0;             % load (1) or not (0) the uncertain vector for best-case problem   
problem_ebro.flag_input.umax_load = 0;             % load (1) or not (0) the uncertain vector for worst-case problem   
problem_ebro.flag_input.fmin_load = 0;             % load (1) or not (0) the value of objective function f(dmin_load, umin_load) 
problem_ebro.flag_input.fmax_load = 0;             % load (1) or not (0) the value of objective function f(dmax_load, umax_load) 
% Type of Output
problem_ebro.flag_output.u_run    = 1;              % Evaluate (1) or not (0) the maximum of f for a fixed design d
problem_ebro.flag_output.dmin_run = 0;              % Evaluate (1) or not (0) the minimum of f for a fixed design u
problem_ebro.flag_output.minmin   = 0;              % Evaluate (1) or not (0) the Best Case Scenario (min-min problem)
problem_ebro.flag_output.minmax   = 0;              % Evaluate (1) or not (0) the Worst Case Scenario (min-max problem)
problem_ebro.flag_output.Belief   = 1;              % Use (1) or not (0) Evidence-Network-Model to reconstruct the Belief curve
problem_ebro.flag_output.Plausibility = 0;          % Use (1) or not (0) Evidence-Network-Model to reconstruct the Plausibility curve
problem_ebro.flag_output.exact_Belief = 0;          % Evaluate (1) or not (0) the exact Belief curve for the fixed design d
problem_ebro.flag_output.exact_Plausibility = 0;    % Evaluate (1) or not (0) the exact Plausibility curve for the fixed design d
problem_ebro.flag_output.plot = 0;                  % plot (1) or not (0) the curves

problem_ebro.objfun     = {@mask_Acta_constr2unconstr_DataVolume};        
problem_ebro.constraint = {[]};
problem_ebro.obj_constr = 0; 


load('min_deterministic');
problem_ebro.input.d_max = d_min_det;
%% ------------------------------------------------------------------------




%% initialise algorithm minmax (META-ALGO)
[ algo_minmax, algo_outer, algo_inner, algo_decomposition ] = init_algorithm(problem_ebro);


algo_decomposition.par.nFeValMax = 1000;
algo_inner.par.nFeValMax = 40000;
%% optimise
[minmin, minmax, LIST, LIST_EXACT] = optimise_minmax_so_decomposition(problem_ebro, algo_inner, algo_outer, algo_minmax, algo_decomposition);



%%  create results directory
if save_flag

        mkdir(savefolder);
        save(strcat(savefolder,'ebro_IAC_max_MASS_d_unconstrained',num2str(algo_decomposition.par.nFeValMax)));
end


