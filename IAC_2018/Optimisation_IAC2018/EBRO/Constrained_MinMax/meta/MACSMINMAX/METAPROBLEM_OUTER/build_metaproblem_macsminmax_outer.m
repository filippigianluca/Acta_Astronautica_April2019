% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [metaproblem] = build_metaproblem_macsminmax_outer(problem_fix_d, local_search)

%% chromosome
dim_d =  problem_fix_d.par_objfun.dim_d;
metaproblem.dim = dim_d;
metaproblem.lb = zeros(1,dim_d);
metaproblem.ub = ones(1,dim_d);


%% function and constraint
metaproblem.par_objfun.problem_fix_d = problem_fix_d;
metaproblem.par_objfun.u_record = [];
metaproblem.par_objfun.local_search = local_search;
metaproblem.par_objfun.objectives = 1:problem_fix_d.n_obj;
% if problem_fix_d.fitnessfcn.obj_constr == 0
metaproblem.objfun = @mask_objfun_macsminmax_outer; % depends on d and par_objfun. Needs to specify par_objfun.(all that is above)
% else
%     metaproblem.objfun = @mask_objfun_constr_macsminmax_outer; 
% end
if isempty(problem_fix_d.par_objfun.constraint{1})
    metaproblem.constraint = [];
else
    metaproblem.constraint = @mask_constraints_macsminmax_outer; %[]; %
end

metaproblem.par_objfun.constraint{1} = problem_fix_d.par_objfun.constraint;




% if isempty(problem_fix_d.par_objfun.constraint{1})
%     metaproblem.par_objfun.mask_constraints = [];
% else
%     metaproblem.par_objfun.mask_constraints = @mask_constraints_macsminmax_outer; %[]; %
% end

% metaproblem.par_objfun.constraints_flag = 0;

%%
% % objective and constraints are defined in different functions
% 
% % Function to optimise
% fitnessfcn.obj       = metaproblem.objfun;
% % Function of constraints
% fitnessfcn.constr    = metaproblem.par_objfun.mask_constraints;
% 
% 
% 
% 
% % Flag to 0: objective and constraints are NOT in the same function
% fitnessfcn.obj_constr = 0;%problem_fix_d.fitnessfcn.obj_constr;
% % How to handle constraints: set to 1 for weighted constraints with fixed
% % weights, or to 0 for penalty with no weights
% fitnessfcn.weighted = problem_fix_d.fitnessfcn.weighted;
% % If the constraints are handled without weights, then define a tolerance
% % for the violation of the equality constraints
% fitnessfcn.ceq_eps = problem_fix_d.fitnessfcn.ceq_eps;
% % Weights for penalty if fitnessfcn.weighted == 1
% fitnessfcn.w_ceq = problem_fix_d.fitnessfcn.w_ceq;
% fitnessfcn.w_c = problem_fix_d.fitnessfcn.w_c;
% 
% 
% metaproblem.fitnessfcn = fitnessfcn;



% Flag to 0: objective and constraints are NOT in the same function
% metaproblem.fitnessfcn.obj_constr = 0;
% Weights for penalty
% metaproblem.fitnessfcn.w_ceq = 100;
% metaproblem.fitnessfcn.w_c = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return