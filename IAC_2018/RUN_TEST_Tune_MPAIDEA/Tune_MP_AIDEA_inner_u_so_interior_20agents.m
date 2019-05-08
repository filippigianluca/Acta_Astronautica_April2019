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
%
%
% tune parameters of MPAIDEA for the inner loop of the minmax problem:
% maximisation over the uncertain vector u of the objective function

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
ALLdir  = CURRENTPATH(1:idcs(end-2)-1);
IAC2018 = CURRENTPATH(1:idcs(end)-1);
if isunix
    Stuff_Danda = strcat(ALLdir,'/utopiae-reliables/Stuff_DANDA');
else
    Stuff_Danda = strcat(ALLdir,'\utopiae-reliables\Stuff_DANDA');
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
% global f
% f = [];

nFeValMax     = 30000;
n_populations = 1;
n_agents      = 20; % prova con 10 e 20
max_LR        = 20;
plots_flag    = 0;
record_vector = 0.001:0.001:1;
fmincon_set = 'interior-point';

nu = 390;

N_fixed_d = 50;

d_fixed = [];

% initialise_problem   = str2func('init_tc_so_unconstr_IAC_2018_tune_MPAIDEA');
initialise_problem   = str2func('init_tc_so_constr_IAC_2018_tune_MPAIDEA_fixed_d');

initialise_algorithm = str2func('init_algo_so_ebro_IAC_2018');





%% ------------------------------------------------------------------------
% RUN
%--------------------------------------------------------------------------
[problem] = initialise_problem();

problem.dim = problem.dim_u;
for k = 1 : size(problem.lb_u{1},1)
    problem.lb(k) = problem.lb_u{1}{k}(1);
    problem.ub(k) = problem.lb_u{1}{k}(end);
end


problem.fitnessfcn.obj = problem.objfun{1};
problem.fitnessfcn.constr = problem.constraints{1};
problem.fitnessfcn.obj_constr = 0;
problem.fitnessfcn.weighted = 0;
problem.fitnessfcn.ceq_eps    = 1e-6;
problem.fitnessfcn.w_ceq = 1000;  % Weights for penalty
problem.fitnessfcn.w_c = 100;

problem.par_objfun.fix.nu = nu;


[~, ~, algo_inner, ~] = initialise_algorithm(problem);

algo_inner.par.nFeValMax     = nFeValMax;
algo_inner.par.n_populations = n_populations;
algo_inner.par.n_agents      = n_agents;
algo_inner.par.max_LR        = max_LR;
algo_inner.par.plots         = plots_flag; 
algo_inner.par.record        = record_vector;  
algo_inner.par.fmincon_set   = fmincon_set;


if isempty(d_fixed)
    % D = rand(N_fixed_d, problem.dim_d);
    load('D_matrix');
else
    D = d_fixed;
end



for i=1:5%size(D,1)
    
    if isempty(d_fixed)
        problem.par_objfun.d = problem.lb_d' + D(i,:).*(problem.ub_d' - problem.lb_d');
    else
        problem.par_objfun.d = d_fixed(i,:);
    end
    
    

    
    % MP-AIDEA optimisation
    [ x, fval, exitflag, output ] = optimise_mpaidea_wrapper(problem,algo_inner.par);
    
save(strcat('tune_constr1obj_inner_30000nfeval_',num2str(n_agents),'agents_D',num2str(i),fmincon_set,'_nu_',num2str(nu)));
    
    
%     nonzeros = [];
%     for k = 1 : size(output.memories_record,1)
%         
%         nonzeros = find(output.memories_record(k,end-1,:)~=0);
%         if ~isempty(nonzeros)
%             [~, b] = min(output.memories_record(k,end-1,nonzeros));
%             memories_record_end(k,:) = output.memories_record(k,:,b);
%         end
% 
%     end
%     
%     recorded_results(:, :, i) = memories_record_end;
% 
%     
% %     N_FF = size(output.population_evolution{1, 1},1);
% %     for j=1:N_FF
% %     FF(i, j) = problem.fitnessfcn.obj(output.population_evolution{1, 1}(j,:), problem.par_objfun);
% %     end
% %     scatter(1:N_FF, FF, 'filled')
% 
% save('Tune_u_IAC2018_2pop_sqp_30000feval')
end
