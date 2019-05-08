% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ----------
%--------------- e-mail: smart@strath.ac.uk -------------------------------
%------------------- Authors: SMART developers team -----------------------
clear all; close all; clc



warning('off','all')
%warning('on','all')

%% ------------------------------------------------------------------------
% add folders to path
%--------------------------------------------------------------------------

repo_folder = '../../';

% addpath(genpath([repo_folder 'smart-o2c/Optimisation']))
% addpath(genpath([repo_folder 'smart-o2c/Other_Tools']))
% addpath(genpath([repo_folder 'smart-o2c/Problems/EBRO']))
rmpath(genpath([repo_folder 'smart-o2c']))
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


%% initialise algorithm minmax (META-ALGO)
[ algo_minmax, algo_outer, algo_inner, algo_decomposition ] = init_algorithm(problem_ebro);

algo_decomposition.par.nFeValMax = 15000;
%% optimise
[minmin, minmax, LIST, LIST_EXACT] = optimise_minmax_so_decomposition(problem_ebro, algo_inner, algo_outer, algo_minmax, algo_decomposition);



%%  create results directory
if save_flag

        mkdir(savefolder);
        save(strcat(savefolder,'ebro_IAC_constr_SO_nfeval_decomposition',num2str(algo_decomposition.par.nFeValMax)));
end


