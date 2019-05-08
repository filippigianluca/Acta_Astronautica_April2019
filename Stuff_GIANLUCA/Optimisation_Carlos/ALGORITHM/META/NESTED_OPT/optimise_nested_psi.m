function [d,fval,exitflag,output] = optimise_nested_psi(problem_minmax, algo_outer, algo_inner, par_minmax)

global nfevalglobal
global history_outer
global surrogate_model

% Rename inputs
n_d = problem_minmax.dim_d;
n_u = problem_minmax.dim_u;
n_obj = problem_minmax.n_obj;
sign_inner = problem_minmax.sign_inner; % 1 for minmax, -1 for minmin

nfevalmax = par_minmax.maxnfeval;

% Build metaproblems
problem_inner = build_metaproblem_nested_inner(problem_minmax);
problem_outer = build_metaproblem_nested_outer_psi(problem_inner, algo_inner);

% Optimise
[duminmax, fval, ~, output_aux] = algo_outer.optimise(problem_outer,algo_outer.par);

nfeval = nfevalglobal; %improve this

frontsize = size(fval,1);
d = duminmax(:,1:n_d).*repmat(problem_minmax.ub_d'-problem_minmax.lb_d',[frontsize,1]) + repmat(problem_minmax.lb_d',[frontsize,1]);
u = cell(1,n_obj);
for obj = 1:n_obj
    map_info = problem_inner.par_objfun.map_u_info{obj};
    for i = 1:frontsize
        if(length(duminmax(i,:))>n_d)
            u{obj}(i,:) = map_affine(duminmax(i,n_d+(obj-1)*n_u+1:n_d+obj*n_u),map_info); %this can be easily vectorized
        else
            u{obj}(i,:) = nan (1,n_u);
        end
    end
end

output.u = u;
output.nfeval = nfeval;
exitflag = 0;

end