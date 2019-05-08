% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [min] = evaluate_min(problem, d, algo_inner)

par_objfun.problem = problem;
par_objfun.objective = problem.n_obj;


n_obj_dec = 1;

par_objfun.problem.flag = -1;  % min
par_objfun.problem.sign = -1;




[meta_algo_inner_max_u] = problem_max_u(problem, n_obj_dec, par_objfun, d, algo_inner);
meta_algo_inner_max_u.par_objfun.flag = -1; %  min

meta_algo_inner_max_u.constraint = @mask_constraint_minimum_fixed_d_so;


meta_algo_inner_max_u.par_objfun.problem_par_objfun{1, 1}.fix = problem.par_objfun{1, 1}.fix;

[ min_u, fmin_u , ~ , ~ ] = meta_algo_inner_max_u.optimise(meta_algo_inner_max_u, meta_algo_inner_max_u.par);

u_true_min = map_affine(min_u, meta_algo_inner_max_u.par_objfun.map_u_info{meta_algo_inner_max_u.par_objfun.objective});




%%%%%%%%%%%%%%%%%
min.f{1} = fmin_u;
min.u{1} = u_true_min;
min.d = d;

end