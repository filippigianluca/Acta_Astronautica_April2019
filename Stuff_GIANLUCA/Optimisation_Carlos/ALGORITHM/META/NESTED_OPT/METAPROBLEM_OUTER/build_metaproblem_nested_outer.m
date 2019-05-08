function [metaproblem] = build_metaproblem_nested_outer(problem_inner,algo_inner)

% chromosome
dim_d =  problem_inner.par_objfun.dim_d;
metaproblem.dim = dim_d;
metaproblem.lb = zeros(1,dim_d+problem_inner.n_obj*problem_inner.dim);
metaproblem.ub = ones(1,dim_d+problem_inner.n_obj*problem_inner.dim);

% function
metaproblem.par_objfun.problem_inner = problem_inner;
metaproblem.par_objfun.objectives = 1:problem_inner.n_obj;
metaproblem.par_objfun.dim_d = dim_d;

metaproblem.par_objfun.algo_inner = algo_inner;
metaproblem.objfun = @mask_objfun_nested_outer; % depends on d and par_objfun.

return
