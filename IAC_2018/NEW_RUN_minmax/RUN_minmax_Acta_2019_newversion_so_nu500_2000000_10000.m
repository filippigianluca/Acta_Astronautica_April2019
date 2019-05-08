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
clear all; close all; clc;


warning('off','all')
repo_folder = '../../';

% addpath(genpath([repo_folder 'smart-o2c/Optimisation']))
% addpath(genpath([repo_folder 'smart-o2c/Other_Tools']))
% addpath(genpath([repo_folder 'smart-o2c/Problems/EBRO']))
addpath(genpath([repo_folder 'smart-o2c/Optimisation']))
addpath(genpath([repo_folder 'utopiae-reliables/IAC_2018']))
addpath(genpath([repo_folder 'utopiae-reliables/Stuff_DANDA']))
rmpath(genpath([repo_folder 'utopiae-reliables/IAC_2018/Optimisation_IAC2018']))


% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);


%% ------------------------------------------------------------------------
% Define Problem
%--------------------------------------------------------------------------


initialise_problem   = str2func('init_tc_so_unconstr_IAC_2018');
initialise_algorithm = str2func('init_algo_so_ebro_IAC_2018');

savefolder = strcat('RESULTS/minmax_so_newversion/');
flag_plot = 1;

%% ------------------------------------------------------------------------
% RUN
%--------------------------------------------------------------------------
[problem_minmax] = initialise_problem();

[algo_minmax, algo_outer, algo_inner, algo_constraint] = initialise_algorithm(problem_minmax);


problem_minmax.par_objfun{1, 1}.fix.nu = 500;
algo_outer.par.nFeValMax  = 10000;
algo_inner.par.nFeValMax  = 10000;
algo_minmax.par_minmax.maxnfeval = 2000000;


% profile on
[ minmin, minmax ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_constraint,algo_minmax.par_minmax);
% [minmin, minmax, LIST, LIST_EXACT] = ebro(problem, algo_inner, algo_outer, algo_minmax, algo_decomposition);
% profile viewer

if flag_plot
    
        mkdir(savefolder);
        save(strcat(savefolder,'Acta_nu500_so_nfeval_minmax_',num2str(algo_minmax.par_minmax.maxnfeval),'nfeval_inner_outer_',num2str(algo_outer.par.nFeValMax)));
end
