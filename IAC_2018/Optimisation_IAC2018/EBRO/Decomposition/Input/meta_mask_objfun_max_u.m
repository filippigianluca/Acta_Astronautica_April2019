% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [varargout] = meta_mask_objfun_max_u(u_to_opt,par_objfun)
% unscales d, u, evaluates the corresponding objfun and multiplies by sign

% To be sure that lb_d, ub_d and d are all row vectors
d = par_objfun.d(:)';
lb_d = par_objfun.lb_d(:)';
ub_d = par_objfun.ub_d(:)';
d_true = lb_d + d.*(ub_d - lb_d);



% d_true = par_objfun.d;                                 % for the moment I don't change d

obj = par_objfun.objective;
map_u_info = par_objfun.map_u_info{obj};

% u_aux = par_objfun.umax;

% u_aux(par_objfun.vars_to_opt)=u_to_opt;
u_true = map_affine(u_to_opt,map_u_info);

%map_u_info = par_objfun.map_u_info{obj};

% u_aux = par_objfun.umax;
% 
% % u_aux(par_objfun.vars_to_opt)=u_to_opt;
% % u_true = map_affine(u_aux,map_u_info);
% 
% %u_true = map_affine(u_to_opt,map_u_info);
% u_aux(par_objfun.vars_to_opt)=u_to_opt;
% u_true = u_aux;


func = par_objfun.objfun{obj};
par_func = par_objfun.problem_par_objfun{1};


F = func(d_true,u_true,par_func);    % -par_objfun.sign*   VERIFICARE 

% masked = -par_objfun.flag*func(d_true,u_true,par_func);    % -par_objfun.sign*   VERIFICARE 


if isstruct(F)
    varargout{1} = -par_objfun.flag*F.f;
    varargout{2} = F.c;
    varargout{3} = F.ceq;
else
    varargout{1} = -par_objfun.flag*F;
end



end