% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------

clear all; close all; clc;


warning('off','all')
%warning('on','all')

%% ------------------------------------------------------------------------
% add folders to path
%--------------------------------------------------------------------------

repo_folder = '../../';

addpath(genpath([repo_folder 'utopiae-reliables/Stuff_DANDA']))
addpath(genpath([repo_folder 'utopiae-reliables/IAC_2018']))
% addpath(genpath([repo_folder 'smart-o2c/Optimisation']))
% rmpath(genpath([repo_folder 'utopiae-reliables/IAC_2018/Optimisation_IAC2018']))

%--------------------------------------------------------------------------



% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);



%% initialisation
init_tc   = str2func(strcat('init_tc_mo_IAC2018'));
init_algo = str2func(strcat('init_algo_mo_IAC2018_mask'));

savefolder = strcat('RESULTS/minmax_mo/');
flag_plot = 1;


disp('MO_')

global nfevalglobal;
nfevalglobal = 0;

%% initialise problem 
[ problem_minmax ] = init_tc();

%% initialise algorithm minmax (META-ALGO).
[ algo_minmax, algo_outer, algo_inner ] = init_algo(problem_minmax);

%% optimise
algo_minmax.optimise = @optimise_mo_ausa_IAC;
[ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);


if flag_plot
    
        mkdir(savefolder);
        save(strcat(savefolder,'IAC_2018_mo_nfeval_minmax_',num2str(algo_minmax.par_minmax.maxnfeval),'nfeval_inner_outer_',num2str(algo_outer.par_au.maxnfeval)));
end
