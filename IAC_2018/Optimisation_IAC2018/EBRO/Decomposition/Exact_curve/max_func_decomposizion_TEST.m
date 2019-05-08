% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [EXACT_FocalElement] = max_func_decomposizion_TEST(in, ii, position_FE, EXACT_FocalElement, algo_inner, problem_inner)

global num_maximization_Belief_exact;
num_maximization_Belief_exact = num_maximization_Belief_exact + 1;


%--------------------------------------------------------------------------
% MAXIMISATION (Belief)
%--------------------------------------------------------------------------

if in.flag_output.exact_Belief
    problem_inner.objfun        = @mask_objfun_max_decomposition;
    problem_inner.constraint    = @mask_constraint_max_decomposition;
    
    % objective and constraints are defined in different functions
    
%     % Function to optimise
%     fitnessfcn.obj       = problem_inner.objfun;
%     % Function of constraints
%     if isempty(problem_inner.par_objfun.constraint{1})
%         fitnessfcn.constraint = [];
%     else
%         fitnessfcn.constraint    = @mask_constraint_max_decomposition;
%     end
    
    
% %--------------------------------------------------------------------------    
%     % Flag to 0: objective and constraints are NOT in the same function
%     fitnessfcn.obj_constr = 0;
%     % How to handle constraints: set to 1 for weighted constraints with fixed
%     % weights, or to 0 for penalty with no weights
%     fitnessfcn.weighted = 0;
%     % If the constraints are handled without weights, then define a tolerance
%     % for the violation of the equality constraints
%     fitnessfcn.ceq_eps = 1e-6;
%     % Weights for penalty if fitnessfcn.weighted == 1
%     fitnessfcn.w_ceq = 100;
%     fitnessfcn.w_c = 100;
%      
%     problem_inner.fitnessfcn = fitnessfcn;
% %--------------------------------------------------------------------------    
    
    problem_inner.par_objfun.d = problem_inner.par_objfun.d_belief;
    problem_inner.par_objfun.u = problem_inner.par_objfun.u_belief;
    problem_inner.par_objfun.lb_d = problem_inner.lb;
    problem_inner.par_objfun.ub_d = problem_inner.ub;
    
    
    problem_inner.par_objfun.maximise = 1;
    
    [ u_max_to_opt_all, fmax_to_opt_all, ~ , output ] = algo_inner.optimise(problem_inner,algo_inner.par);
    
    

    [fmax_to_opt, N_fmax_to_opt] = min(fmax_to_opt_all);
    u_max_to_opt = u_max_to_opt_all(N_fmax_to_opt, :);
    
    
    EXACT_FocalElement{1,ii}.upper_u = u_max_to_opt;
    EXACT_FocalElement{1,ii}.upper_f = -fmax_to_opt;
    
    EXACT_FocalElement{1,ii}.output_mpaidea = output;
end



%--------------------------------------------------------------------------
% MINIMISATION (Plausibility)
%--------------------------------------------------------------------------
if in.flag_output.exact_Plausibility
    
    problem_inner.objfun = @mask_objfun_min_decomposition;
    
    
    % objective and constraints are defined in different functions
    
    % Function to optimise
    fitnessfcn.obj       = problem_inner.objfun;
    % Function of constraints
    if isempty(problem_inner.par_objfun.constraint{1})
        fitnessfcn.constr = [];
    else
        fitnessfcn.constr    = @mask_constraint_max_decomposition;
    end
    
    
% %--------------------------------------------------------------------------    
%     % Flag to 0: objective and constraints are NOT in the same function
%     fitnessfcn.obj_constr = 0;
%     % How to handle constraints: set to 1 for weighted constraints with fixed
%     % weights, or to 0 for penalty with no weights
%     fitnessfcn.weighted = 0;
%     % If the constraints are handled without weights, then define a tolerance
%     % for the violation of the equality constraints
%     fitnessfcn.ceq_eps = 1e-6;
%     % Weights for penalty if fitnessfcn.weighted == 1
%     fitnessfcn.w_ceq = 100;
%     fitnessfcn.w_c = 100;
%     
%     problem_inner.fitnessfcn = fitnessfcn;
% %--------------------------------------------------------------------------

    problem_inner.par_objfun.d = problem_inner.par_objfun.d_plausibility;
    problem_inner.par_objfun.u = problem_inner.par_objfun.u_plausibility;
    problem_inner.par_objfun.lb_d = problem_inner.lb;
    problem_inner.par_objfun.ub_d = problem_inner.ub;
    
    problem_inner.par_objfun.maximise = 0;
    
    [ u_min_to_opt_all, fmin_to_opt_all , ~ , ~ ] = algo_inner.optimise(problem_inner,algo_inner.par);

    
    [fmin_to_opt, N_fmin_to_opt] = min(fmin_to_opt_all);
    u_min_to_opt = u_min_to_opt_all(N_fmin_to_opt, :);

    
    EXACT_FocalElement{1,ii}.downer_u = u_min_to_opt;
    EXACT_FocalElement{1,ii}.downer_f = fmin_to_opt;
end



bpa = 1;
for k = 1:problem_inner.dim
    bpa = bpa*problem_inner.bpa{k,1}(position_FE(k));  % problem_inner.bpa{1,1}{k,1}(position_FE(k));
end

EXACT_FocalElement{1,ii}.bpa = bpa;
EXACT_FocalElement{1,ii}.n_FE = ii;
EXACT_FocalElement{1,ii}.lb = problem_inner.lb;
EXACT_FocalElement{1,ii}.ub = problem_inner.ub;


end