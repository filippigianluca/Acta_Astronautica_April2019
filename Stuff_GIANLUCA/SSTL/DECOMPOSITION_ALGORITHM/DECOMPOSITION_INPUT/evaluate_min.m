function [min] = evaluate_min(problem, d, algo_inner)

par_objfun.problem = problem;
par_objfun.objective = problem.n_obj;


n_obj_dec = 1;

par_objfun.problem.flag = -1;  % min

[meta_algo_inner_max_u] = problem_max_u(problem, n_obj_dec, par_objfun, d, algo_inner);
meta_algo_inner_max_u.par_objfun.flag = -1; %  min
[ min_u, fmin_u , ~ , ~ ] = meta_algo_inner_max_u.optimise(meta_algo_inner_max_u, meta_algo_inner_max_u.par);

u_true_min = map_affine(min_u, meta_algo_inner_max_u.par_objfun.map_u_info{meta_algo_inner_max_u.par_objfun.objective});


%%%%%%%%%%%%%%%%%
min.f = fmin_u;
min.u{1} = u_true_min;
min.d = d;

end