% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [algo_minmax, algo_outer, algo_inner] = init_algo_mo_IAC2018_mask(minmax_problem)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meta-algorithm
nfevalmax_total = 5e6;
local_search_validation = false;
% outer loop
nfevalmax_outer        = 30000;% 100*minmax_problem.dim_d;
population_outer       = 10 + minmax_problem.dim_d;
outer_front_size       = 7;   % integer [2,inf]
ratio_surrogate_assist = 0;  % real [0,1]
local_search_outer     = false;
% inner loop
nfevalmax_inner      = [20000 20000]; %[600*minmax_problem.dim_u, 400*minmax_problem.dim_u];
population_inner     = max(4,minmax_problem.dim_u);
local_search_inner   = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OUTPUT
minmax_problem.maxnfeval = nfevalmax_total;
[algo_minmax, algo_outer, algo_inner] = init_algo_mo_IAC2018(minmax_problem);

algo_minmax.par_minmax.maxnfeval = nfevalmax_total;
algo_minmax.par_minmax.local_search_flags.validation = local_search_validation;
algo_minmax.par_minmax.local_search_flags.outer = local_search_outer;

if (size(nfevalmax_inner) == [1 1])
    nfevalmax_inner = repmat(nfevalmax_inner,[1,minmax_problem.n_obj]);
elseif (size(nfevalmax_inner) ~= [1,minmax_problem.n_obj])
    error('nfevalmax_inner wrong format')
end

if (size(population_inner) == [1 1])
    population_inner = repmat(population_inner,[1,minmax_problem.n_obj]);
elseif (size(population_inner) ~= [1,minmax_problem.n_obj])
    error('population_inner wrong format')
end

for obj = 1:minmax_problem.n_obj
    algo_inner{obj}.par.nFeValMax = nfevalmax_inner(1,obj);
    algo_inner{obj}.par.n_agents = population_inner(1,obj);

end


if (size(local_search_inner) == [1 1])
    local_search_inner = repmat(local_search_inner,[1,minmax_problem.n_obj]);
elseif (size(local_search_inner) ~= [1,minmax_problem.n_obj])
    error('local_search_inner wrong format')
end

algo_minmax.par_minmax.local_search_flags.inner = local_search_inner;

algo_outer.par_au.maxnfeval = nfevalmax_outer;
algo_outer.par_au.popsize = population_outer;

algo_outer.par_sa.max_arch_out = round(ratio_surrogate_assist*outer_front_size);
algo_outer.par_au.max_arch_out = outer_front_size - algo_outer.par_sa.max_arch_out;

if(algo_outer.par_au.max_arch_out==1)
    algo_outer.par_au.max_arch_out = 0;
    algo_outer.par_sa.max_arch_out = algo_outer.par_au.max_arch_out + 1;
    warning(strcat('increased the ratio of surrogate assist to ',' ',...
        num2str(algo_outer.par_sa.max_arch_out/(algo_outer.par_sa.max_arch_out+algo_outer.par_au.max_arch_out)),...
        ' to avoid having just one individual in the minimisation front'))
elseif(algo_outer.par_sa.max_arch_out==1)
    algo_outer.par_sa.max_arch_out = 2;
    algo_outer.par_au.max_arch_out = algo_outer.par_au.max_arch_out - 1;
    warning(strcat('increased the ratio of surrogate assist to ',' ',...
        num2str(algo_outer.par_sa.max_arch_out/(algo_outer.par_sa.max_arch_out+algo_outer.par_au.max_arch_out)),...
        ' to avoid having just one individual in the surrogate front'))
end


if(algo_outer.par_au.maxnfeval == 0 || algo_outer.par_au.max_arch_out == 0)
    algo_outer.par_au.use = false;
end

if(algo_outer.par_sa.max_arch_out == 0)
    algo_outer.par_sa.use = false;
end

if (~algo_outer.par_au.use && ~algo_outer.par_sa.use)
    error('no minimisation at all in outer loop; you probably set bot nfevalmax_outer and ratio_surrogate_assist to zero')
end



algo_outer.par_au.func  = minmax_problem.objfun;
algo_outer.par_au.cfunc = minmax_problem.constraints;

algo_outer.par_sa.func  = minmax_problem.objfun;
algo_outer.par_sa.cfunc = minmax_problem.constraints;


return