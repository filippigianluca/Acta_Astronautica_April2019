function [masked, du] = mask_objfun_nested_outer(d,par_objfun)


par_objfun.problem_inner.par_objfun.d = d; % tell metaproblem inner what d to fix
masked = [];
u = [];
for obj = par_objfun.objectives
    par_objfun.problem_inner.par_objfun.objective = obj;                                                                                       % tell metaproblem what objective to optimise
    [ u_inner, f_inner , ~ , output_aux ] = par_objfun.algo_inner{obj}.optimise(par_objfun.problem_inner,par_objfun.algo_inner{obj}.par);      % optimise
    u = [u, u_inner];
    masked = [masked, -par_objfun.problem_inner.par_objfun.sign*f_inner];
end

du = [d(1:par_objfun.dim_d) u];
return
