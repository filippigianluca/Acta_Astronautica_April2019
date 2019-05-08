function [metaproblem] = build_metaproblem_nested_outer_psi(problem_inner,algo_inner)

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

% surrogate
% These are hardcoded (HC) at the moment but should really come from par_minmax for flexibility
    surrogate.method = 'kriging';
    surrogate.corrfun = @corrgauss;
    surrogate.regrfun = @regpoly0;
    surrogate.training = str2func([lower(surrogate.method) '_training']);
    surrogate.predictor = str2func([lower(surrogate.method) '_predictor']);
    % surrogate.model = []; % the model is a global variable for now

    surrogate.update_frequency = 10;
    surrogate.max_set_size = 100;
    surrogate.set_size_minima = 15; % when building the surrogate it will try to pick this subset of the dataset amongst the best solutions so far.

metaproblem.par_objfun.surrogate = surrogate;

metaproblem.objfun = @mask_objfun_nested_outer_psi; % depends on d and par_objfun.

return
