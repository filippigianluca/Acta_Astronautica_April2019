function [metaproblem] = build_metaproblem_ideaminmax_s_inner(problem)

% chromosome
dim_u =  problem.dim_u;
metaproblem.dim = dim_u;
metaproblem.lb = zeros(1,dim_u);
metaproblem.ub = ones(1,dim_u);

metaproblem.par_objfun.surrogate = [];

metaproblem.par_objfun.sign = problem.sign_inner; %dunno if necessary
metaproblem.par_objfun.ymin = [];

metaproblem.objfun = @mask_objfun_ideaminmax_s_inner; %depends on u and par_objfun. Needs to specify par_objfun.surrogate.model, par_objfun.ymin

return