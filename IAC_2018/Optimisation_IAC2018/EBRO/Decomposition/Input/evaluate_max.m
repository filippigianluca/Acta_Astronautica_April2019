% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [max] = evaluate_max(problem, d, algo_inner)


% normalise fixed design d
d = d(:)'; lb_d = problem.lb_d(:)'; ub_d = problem.ub_d(:)'; 
d = (d-lb_d)./(ub_d - lb_d);

par_objfun.problem = problem;
par_objfun.objective = problem.n_obj;

n_obj_dec = 1;

[meta_algo_inner_max_u] = problem_max_u(problem, n_obj_dec, par_objfun, d, algo_inner);
meta_algo_inner_max_u.par_objfun.flag = 1; 





meta_algo_inner_max_u.par_objfun.problem_par_objfun{n_obj_dec}.fix = problem.par_objfun{n_obj_dec}.fix;


[ max_u_pop, fmax_u_pop, ~ , ~ ] = meta_algo_inner_max_u.optimise(meta_algo_inner_max_u, meta_algo_inner_max_u.par);


[fmax_u, N_fmax_u] = min(fmax_u_pop);
max_u = max_u_pop(N_fmax_u, :);


u_true = map_affine(max_u, meta_algo_inner_max_u.par_objfun.map_u_info{meta_algo_inner_max_u.par_objfun.objective});


max.u{1} = u_true;
max.f{1} = -fmax_u;
max.d = lb_d + d.*(ub_d-lb_d);

end