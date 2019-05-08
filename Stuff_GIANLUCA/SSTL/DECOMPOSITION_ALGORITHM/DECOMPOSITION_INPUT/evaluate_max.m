function [max] = evaluate_max(problem, d, algo_inner)

par_objfun.problem = problem;
par_objfun.objective = problem.n_obj;

n_obj_dec = 1;

[meta_algo_inner_max_u] = problem_max_u(problem, n_obj_dec, par_objfun, d, algo_inner);
meta_algo_inner_max_u.par_objfun.flag = 1; %  max
[ max_u, fmax_u , ~ , ~ ] = meta_algo_inner_max_u.optimise(meta_algo_inner_max_u, meta_algo_inner_max_u.par);

u_true = map_affine(max_u, meta_algo_inner_max_u.par_objfun.map_u_info{meta_algo_inner_max_u.par_objfun.objective});

%%%%%%%%%%%%%%%%%%%%
max.u{1} = u_true;
max.f = -fmax_u;
max.d = d;

end