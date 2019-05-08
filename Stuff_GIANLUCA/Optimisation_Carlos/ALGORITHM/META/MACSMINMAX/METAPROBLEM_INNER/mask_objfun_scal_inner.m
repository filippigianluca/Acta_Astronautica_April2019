function [masked] = mask_objfun_scal_inner(u,par_objfun)
% unscales d, u, evaluates the corresponding objfun and multiplies by sign
% this script needs to change a lot if interfaced to a vectorial function

d_true = par_objfun.lb_d' + par_objfun.d.*(par_objfun.ub_d' - par_objfun.lb_d');
% obj = par_objfun.objective;
n_obj = length(par_objfun.lambdas);


% in a vectorial function all below here changes!
f = nan(1,n_obj);
for obj = 1:n_obj
    if par_objfun.lambdas(1,obj) >= 1e-8
        u_true = map_affine(u,par_objfun.map_u_info{obj});
        % history_aux = [history_aux, u_true];
        func = par_objfun.objfun{obj};
        par_func = par_objfun.problem_par_objfun{obj};
        f(1,obj) = func(d_true,u_true,par_func);
    else
        f(1,obj) = nan;
    end 
end

% history_aux = [history_aux, f];

% global history; history = [history; history_aux];
% global history; history = [history; d u f];


f_scal = -par_objfun.sign*(par_objfun.lambdas.*(f-par_objfun.utopia)./par_objfun.frontspans);
masked = min(f_scal);

return