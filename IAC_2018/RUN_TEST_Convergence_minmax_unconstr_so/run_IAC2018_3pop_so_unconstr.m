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
%warning('on','all')
%--------------------------------------------------------------------------
% add folders to path
%--------------------------------------------------------------------------
CURRENTPATH=pwd;
if isunix
    idcs   = strfind(CURRENTPATH,'/');
else
    idcs   = strfind(CURRENTPATH,'\');
end
ALLdir = CURRENTPATH(1:idcs(end-2)-1);
IAC2018 = CURRENTPATH(1:idcs(end)-1);
if isunix
    Stuff_Danda = strcat(ALLdir,'/UTOPIAE_RELIABLES/Stuff_DANDA');
else
    Stuff_Danda = strcat(ALLdir,'\UTOPIAE_RELIABLES\Stuff_DANDA');
end
rmpath(genpath(ALLdir));
addpath(genpath(IAC2018));
addpath(genpath(Stuff_Danda));

%--------------------------------------------------------------------------
% Reset random numbers generator
%--------------------------------------------------------------------------
seed = 1;
s = RandStream('mt19937ar','Seed', seed);
RandStream.setGlobalStream(s);


%% ------------------------------------------------------------------------
% Define Problem
%--------------------------------------------------------------------------
initialise_problem   = str2func('init_tc_so_unconstr_IAC_2018');
initialise_algorithm = str2func('init_algo_so_ebro_IAC_2018');


N_minmax      = 1e6;
N_InnerOuter  = 20000;
n_populations = 3;
max_LR        = [];
plots_flag    = 0;
record_vector = 1;
fmincon_set = 'sqp';




%% ------------------------------------------------------------------------
% RUN
%--------------------------------------------------------------------------
[problem] = initialise_problem();

[algo_minmax, algo_outer, algo_inner, algo_decomposition] = initialise_algorithm(problem);



%==========================================================================
algo_minmax.par_minmax.maxnfeval = N_minmax;

algo_inner.par.nFeValMax = N_InnerOuter;
algo_inner.par.n_populations = n_populations;
algo_inner.par.max_LR        = max_LR;
algo_inner.par.plots         = plots_flag;
algo_inner.par.record        = record_vector;
algo_inner.par.fmincon_set   = fmincon_set;

algo_outer.par.nFeValMax = N_InnerOuter;
algo_outer.par.n_populations = n_populations;
algo_outer.par.max_LR        = max_LR;
algo_outer.par.plots         = plots_flag;
algo_outer.par.record        = record_vector;
algo_outer.par.fmincon_set   = fmincon_set;
%==========================================================================



%% RUN
[minmin, minmax, LIST, LIST_EXACT] = ebro(problem, algo_inner, algo_outer, algo_minmax, algo_decomposition);


%% save
save(strcat('IAC2018_nfeval_',num2str(N_minmax),'nInnerOuter_',num2str(N_InnerOuter),'npop_',num2str(n_populations),fmincon_set));

