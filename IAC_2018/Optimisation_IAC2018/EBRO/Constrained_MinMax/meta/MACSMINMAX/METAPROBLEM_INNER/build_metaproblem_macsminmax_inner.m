% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [metaproblem] = build_metaproblem_macsminmax_inner(problem)

% for objective function mask
metaproblem.par_objfun.dim_d = problem.dim_d;
metaproblem.par_objfun.lb_d = problem.lb_d;
metaproblem.par_objfun.ub_d = problem.ub_d;
metaproblem.n_obj = problem.n_obj;


constraint_flag = 0;
for obj = 1:problem.n_obj
    metaproblem.par_objfun.objfun{obj}             = problem.objfun{obj};
    metaproblem.par_objfun.constraint{obj}         = problem.constraint{obj};
    metaproblem.par_objfun.obj_constr              = problem.obj_constr;
    metaproblem.par_objfun.problem_par_objfun{obj} = problem.par_objfun{obj};
    % metaproblem.par_objfun.lb_u{obj} = problem.lb_u{obj};
    % metaproblem.par_objfun.ub_u{obj} = problem.ub_u{obj};
    metaproblem.par_objfun.map_u_info{obj}         = get_map_info(problem.lb_u{obj}, problem.ub_u{obj});
    if ~isempty(problem.constraint{obj})
        constraint_flag = 1;
    end
end



%% chromosome
dim_u =  problem.dim_u;
metaproblem.dim = dim_u;
metaproblem.lb = zeros(1,dim_u);
metaproblem.ub = ones(1,dim_u);


%% function and constraint
metaproblem.par_objfun.sign = problem.sign_inner;
% if problem.func_constraints == 0
metaproblem.objfun = @mask_objfun_macsminmax_inner; % depends on u and par_objfun. Needs to specify par_objfun.d and par_objfun.objective
% else
%     metaproblem.objfun = @mask_objfun_constraint_macsminmax_inner; 
% end

if constraint_flag == 0
    metaproblem.constraint = [];
else
    metaproblem.constraint = @mask_constraint_macsminmax_inner;  
end


return