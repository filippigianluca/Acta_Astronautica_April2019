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
ALLdir    = CURRENTPATH(1:idcs(end-2)-1);
MY_GITHUB = CURRENTPATH(1:idcs(end-1)-1);
if isunix
    Stuff_Danda = strcat(MY_GITHUB,'/UTOPIAE_RELIABLES/Stuff_DANDA');
    IAC2018     = strcat(MY_GITHUB,'/UTOPIAE_RELIABLES/IAC_2018');
else
    Stuff_Danda = strcat(MY_GITHUB,'\UTOPIAE_RELIABLES\Stuff_DANDA');
    IAC2018     = strcat(MY_GITHUB,'\UTOPIAE_RELIABLES\IAC_2018');
end

rmpath(genpath(ALLdir));
addpath(genpath(IAC2018));
addpath(genpath(Stuff_Danda));


% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);


%% ------------------------------------------------------------------------
% Define Problem
%--------------------------------------------------------------------------


initialise_problem   = str2func('init_tc_so_constr_IAC_2018');
initialise_algorithm = str2func('init_algo_so_constr_ebro_IAC_2018');

savefolder = strcat('RESULTS/minmax_so_constr/');
flag_plot = 1;

%% ------------------------------------------------------------------------
% RUN
%--------------------------------------------------------------------------
[problem] = initialise_problem();

[algo_minmax, algo_outer, algo_inner, algo_decomposition] = initialise_algorithm(problem);

%%%
algo_minmax.optimise = @optimise_so_run_constr_nu400_interior;
problem.par_objfun{1, 1}.fix.nu = 400;
problem.fix.nu = 400;
algo_minmax.par_minmax.maxnfeval = 1000000;
algo_inner.par.nFeValMax = 10000;
algo_outer.par.nFeValMax = 10000;
algo_inner.par.fmincon_set = 'interior-point';
algo_outer.par.fmincon_set = 'interior-point';
%%%

[minmin, minmax, LIST, LIST_EXACT] = ebro(problem, algo_inner, algo_outer, algo_minmax, algo_decomposition);


if flag_plot
    
        mkdir(savefolder);
        save(strcat(savefolder,'nu400_interior_IAC_2018_so_constr_nfeval_minmax_',num2str(algo_minmax.par_minmax.maxnfeval),'nfeval_inner_outer_',num2str(algo_outer.par.nFeValMax)));
end
