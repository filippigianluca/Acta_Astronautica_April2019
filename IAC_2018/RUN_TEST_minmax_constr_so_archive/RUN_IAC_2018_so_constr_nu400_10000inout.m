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



%% ------------------------------------------------------------------------
% add folders to path
%--------------------------------------------------------------------------

repo_folder = '../../../';

addpath(genpath([repo_folder 'utopiae-reliables/Stuff_DANDA']))
addpath(genpath([repo_folder 'utopiae-reliables/IAC_2018']))
addpath(genpath([repo_folder 'smart-o2c/Optimisation']))
rmpath(genpath([repo_folder 'utopiae-reliables/IAC_2018/Optimisation_IAC2018']))


%--------------------------------------------------------------------------
% warning('off','all')
% %warning('on','all')
% %--------------------------------------------------------------------------
% % add folders to path
% %--------------------------------------------------------------------------
% CURRENTPATH=pwd;
% if isunix
%     idcs   = strfind(CURRENTPATH,'/');
% else
%     idcs   = strfind(CURRENTPATH,'\');
% end
% ALLdir    = CURRENTPATH(1:idcs(end-3)-1);
% MY_GITHUB = CURRENTPATH(1:idcs(end-2)-1);
% if isunix
%     Stuff_Danda = strcat(MY_GITHUB,'/utopiae-reliables/Stuff_DANDA');
%     IAC2018     = strcat(MY_GITHUB,'/utopiae-reliables/IAC_2018');
% else
%     Stuff_Danda = strcat(MY_GITHUB,'\utopiae-reliables\Stuff_DANDA');
%     IAC2018     = strcat(MY_GITHUB,'\utopiae-reliables\IAC_2018');
% end
% 
% rmpath(genpath(ALLdir));
% addpath(genpath(IAC2018));
% addpath(genpath(Stuff_Danda));


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

[algo_minmax, algo_outer, algo_inner] = initialise_algorithm(problem);

%%%
algo_minmax.optimise = @optimise_so_run_constr_nu400_10000inout;
problem.par_objfun{1, 1}.fix.nu = 400;
problem.fix.nu = 400;
algo_minmax.par_minmax.maxnfeval = 5000000;
algo_inner.par.nFeValMax = 10000;
algo_outer.par.nFeValMax = 10000;
algo_inner.par.fmincon_set = 'interior-point';
algo_outer.par.fmincon_set = 'interior-point';

% load('nu390_IAC_2018_so_constr_nfeval_minmax_2000000nfeval_inner_outer_20000','minmax');
algo_minmax.par_minmax.Ad = [];% cell2mat(minmax.output.ARCHIVE{1, 1}(2:end,1));
algo_minmax.par_minmax.Au = [];%  cell2mat(minmax.output.ARCHIVE{1, 1}(2:end,2));
algo_minmax.par_minmax.Ac = [];%  cell2mat(minmax.output.ARCHIVE{1, 1}(2:end,3));
algo_minmax.par_minmax.Af = [];%  cell2mat(minmax.output.ARCHIVE{1, 1}(2:end,4));
algo_minmax.par_minmax.C = [];%  cell2mat(minmax.output.ARCHIVE{1, 1}(2:end,5));
algo_minmax.par_minmax.f = [];%  cell2mat(minmax.output.ARCHIVE{1, 1}(2:end,6));
algo_minmax.par_minmax.N = [];%  cell2mat(minmax.output.ARCHIVE{1, 1}(2:end,7));
%%%

[ minmin, minmax ] = algo_minmax.optimise(problem,algo_outer,algo_inner,algo_minmax.par_minmax);
% [minmin, minmax, LIST, LIST_EXACT] = ebro(problem, algo_inner, algo_outer, algo_minmax, algo_decomposition);


if flag_plot
    
        mkdir(savefolder);
        save(strcat(savefolder,'nu400_IAC_2018_so_constr_nfeval_minmax_',num2str(algo_minmax.par_minmax.maxnfeval),'nfeval_inner_outer_',num2str(algo_outer.par.nFeValMax)));
end
