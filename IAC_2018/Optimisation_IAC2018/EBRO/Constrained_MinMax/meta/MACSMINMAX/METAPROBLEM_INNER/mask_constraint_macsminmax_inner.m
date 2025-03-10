% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [masked, ceq] = mask_constraint_macsminmax_inner(u,par_objfun)
% unscales d, u, evaluates the corresponding objfun and multiplies by sign

% To be sure that lb_d, ub_d and d are all row vectors
d = par_objfun.d(:)';
lb_d = par_objfun.lb_d(:)';
ub_d = par_objfun.ub_d(:)';
d_true = lb_d + d.*(ub_d - lb_d);

if par_objfun.sign == -1
    par_objfun.d = u(1:par_objfun.dim_d);
    u = u(1+par_objfun.dim_d:end);
    par_objfun.objective = 1;
end

%     d_true = par_objfun.lb_d' + par_objfun.d.*(par_objfun.ub_d' - par_objfun.lb_d');
    obj = par_objfun.objective;
    map_u_info = par_objfun.map_u_info{obj};

    u_true = map_affine(u,map_u_info);


    func = par_objfun.constraint{obj};
    par_func = par_objfun.problem_par_objfun{obj};

    masked = par_objfun.sign*func(d_true,u_true,par_func);


ceq = [];


return