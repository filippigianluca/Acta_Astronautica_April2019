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
ALLdir  = CURRENTPATH(1:idcs(end-4)-1);
IAC2018 = CURRENTPATH(1:idcs(end-2)-1);
if isunix
    IAC2018_END = strcat(IAC2018,'/FORMULATION_END');
    Stuff_Danda = strcat(ALLdir,'/UTOPIAE_RELIABLES/Stuff_DANDA');
else
    IAC2018_END = strcat(IAC2018,'\FORMULATION_END');
    Stuff_Danda = strcat(ALLdir,'\UTOPIAE_RELIABLES\Stuff_DANDA');
end

rmpath(genpath(ALLdir));
addpath(genpath(IAC2018_END));
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
% global f
% f = [];

nFeValMax     = 30000;
n_populations = 2;
% n_agents      = 12;
max_LR        = [];
plots_flag    = 1;
record_vector = 0.2:0.1:1;
fmincon_set = 'sqp';

N_fixed_u = 50;

initialise_problem   = str2func('init_tc_so_unconstr_IAC_2018_tune_fixed_d_or_u');
initialise_algorithm = str2func('init_algo_so_ebro_IAC_2018');





%% ------------------------------------------------------------------------
% RUN
%--------------------------------------------------------------------------
[problem] = initialise_problem();

problem.dim = problem.dim_d;
problem.lb = problem.lb_d';
problem.ub = problem.ub_d';

problem.objfun =     {@IAC_function_tune_fix_u};

problem.fitnessfcn.obj = problem.objfun{1};
problem.fitnessfcn.constr = [] ;
problem.fitnessfcn.obj_constr = 0;
problem.fitnessfcn.weighted = 0;


[~, ~, algo_inner, ~] = initialise_algorithm(problem);

algo_inner.par.nFeValMax     = nFeValMax;
algo_inner.par.n_populations = n_populations;
% algo_inner.par.n_agents      = n_agents;
algo_inner.par.max_LR        = max_LR;
algo_inner.par.plots         = plots_flag; 
algo_inner.par.record        = record_vector;  
algo_inner.par.fmincon_set   = fmincon_set;


% U = rand(N_fixed_u, problem.dim_u);
load('U_matrix');

for i=1:N_fixed_u

    map_info = get_map_info(problem.lb_u{1}, problem.ub_u{1});
    problem.par_objfun.u = map_affine(U(i,:), map_info);
    
    
    % MP-AIDEA optimisation
    [ x, fval, exitflag, output ] = optimise_mpaidea_wrapper(problem,algo_inner.par);
   
    nonzeros = [];
    for k = 1 : size(output.memories_record,1)
        
        nonzeros = find(output.memories_record(k,end-1,:)~=0);
        if ~isempty(nonzeros)
            [~, b] = min(output.memories_record(k,end-1,nonzeros));
            memories_record_end(k,:) = output.memories_record(k,:,b);
        end

    end
    
    recorded_results(:, :, i) = memories_record_end;
    
%     N_FF = size(output.population_evolution{1, 1},1);
%     for j=1:N_FF
%     FF(i, j) = problem.fitnessfcn.obj(output.population_evolution{1, 1}(j,:), problem.par_objfun);
%     end
%     scatter(1:N_FF, FF, 'filled')
save('Tune_IAC2018_2pop_sqp_30000feval')

end



