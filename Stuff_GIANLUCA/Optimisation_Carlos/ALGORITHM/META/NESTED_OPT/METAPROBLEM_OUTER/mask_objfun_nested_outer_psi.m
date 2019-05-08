function [masked, du] = mask_objfun_nested_outer_psi(d,par_objfun)


par_objfun.problem_inner.par_objfun.d = d; % tell metaproblem inner what d to fix
masked = [];
u = [];
u0 = [];
global surrogate_model;
global history_outer;

if (~isempty(surrogate_model))
	[u0,~] = par_objfun.surrogate.predictor(d, surrogate_model);
end

for obj = par_objfun.objectives
    par_objfun.problem_inner.par_objfun.objective = obj;                                                                                       % tell metaproblem what objective to optimise
    if (~isempty(u0))
    	par_objfun.algo_inner{obj}.par.initial_population=u0((obj-1)*par_objfun.problem_inner.dim+1 : obj*par_objfun.problem_inner.dim);
    end
    [ u_inner, f_inner , ~ , output_aux ] = par_objfun.algo_inner{obj}.optimise(par_objfun.problem_inner,par_objfun.algo_inner{obj}.par);      % optimise
    u = [u, u_inner];
    masked = [masked, -par_objfun.problem_inner.par_objfun.sign*f_inner];
end

du = [d(1:par_objfun.dim_d) u];


history_outer = [history_outer; du masked];

global nfevalouter
nfevalouter = nfevalouter + 1;

n_obj = max(par_objfun.objectives);

if (mod(nfevalouter, par_objfun.surrogate.update_frequency) == 0)

	dataset = [];

	dom = dominance(history_outer(:,end-n_obj+1:end),0);
	[dom,idx_sort] = sort(dom);
	history_outer = history_outer(idx_sort,:); %the history is now sorted by dominance

	[~,idx_unique] = unique(round(1e8*history_outer(:,1:par_objfun.dim_d)),'rows');
    history_outer = history_outer(idx_unique,:);
    dom = dom(idx_unique); % the history is now cleant of "repeated individuals" and the dominance index is corresponding

    size_history = size(history_outer,1);
    available = logical(ones(1,size_history));

    if (size_history < par_objfun.surrogate.max_set_size)
    	dataset = history_outer;
    else
    	
	    idx_front = dom == 0;

		if (sum(idx_front) > par_objfun.surrogate.set_size_minima)
			error('NEED TO IMPLEMENT AN ARCHIVE SHRINK HERE'); % will never happen in SO
		else
			take = 1:par_objfun.surrogate.set_size_minima;
			available(take) = 0;
			dataset = [dataset; history_outer(take,:)];
		end

		dataset = [dataset; datasample(history_outer(available,:),par_objfun.surrogate.max_set_size - par_objfun.surrogate.set_size_minima,1,'Replace',false)];
	end

	umax = dataset(:,par_objfun.dim_d+1:end-n_obj);
	dset = dataset(:,1:par_objfun.dim_d);

	surrogate_model = par_objfun.surrogate.training(dset,umax,par_objfun.surrogate);

end

	

return
