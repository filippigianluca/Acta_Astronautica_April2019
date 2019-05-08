function [metaproblem] = build_metaproblem_ei_outer(problem)

% chromosome
dim_d =  problem.dim_d;
metaproblem.dim = dim_d;
metaproblem.lb = zeros(1,dim_d);
metaproblem.ub = ones(1,dim_d);

% surrogate
% These are hardcoded (HC) at the moment but should really come from par_minmax for flexibility
    surrogate.method = 'kriging';
    surrogate.corrfun = @corrgauss;
    surrogate.regrfun = @regpoly0;
    surrogate.training = str2func([lower(surrogate.method) '_training']);
    surrogate.predictor = str2func([lower(surrogate.method) '_predictor']);
    surrogate.model = [];
    surrogate.max_set_size = 150;
    surrogate.set_size_minima = 15; % when building the surrogate it will try to pick this subset of the dataset amongst the best solutions so far.

%EI
	if problem.n_obj == 1
	        surrogate.indicator = str2func([lower(surrogate.method) '_EI']);
    elseif problem.n_obj == 2
        surrogate.indicator = str2func([lower(surrogate.method) '_EIaug']);
    else
        error('We do not have kriging indicators for many-objective (yet)' );
    end
	

metaproblem.par_objfun.ymin = [];
metaproblem.par_objfun.surrogate = surrogate;
metaproblem.par_objfun.n_obj = problem.n_obj;


% function
metaproblem.objfun = @mask_objfun_ei_outer; % depends on d and par_objfun. Needs to specify par_objfun.surrogate.model
return