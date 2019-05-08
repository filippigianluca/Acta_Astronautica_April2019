% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [decomposition] = max_min_func_decomposition(i, ii, in, position_FE, decomposition, algo_inner, problem_decomposition)


bpa = 1;
for k = 1:problem_decomposition.dim
    bpa = bpa*problem_decomposition.bpa{k,1}(position_FE(k));  % product of the bpa of the coupled components
end

decomposition{1,i-in.num_functions}.FocalElement{1,ii}.bpa = bpa;
decomposition{1,i-in.num_functions}.FocalElement{1,ii}.n_FE = ii;
decomposition{1,i-in.num_functions}.FocalElement{1,ii}.lb = problem_decomposition.lb;
decomposition{1,i-in.num_functions}.FocalElement{1,ii}.ub = problem_decomposition.ub;
decomposition{1,i-in.num_functions}.FocalElement{1,ii}.position = position_FE;

%--------------------------------------------------------------------------
% MAX
%--------------------------------------------------------------------------
if in.flag_output.Belief   %in.output == 0 || in.output == 2
    
    
    
    problem_decomposition.objfun = @mask_objfun_max_decomposition;
    
    % objective and constraints are defined in different functions
    % Function to optimise
%     fitnessfcn.obj       = problem_inner.objfun;
    % Function of constraints
    if isempty(problem_decomposition.par_objfun.constraint{1})
        problem_decomposition.constraint    = [];
    else
        problem_decomposition.constraint    = @mask_constraint_max_decomposition;
    end
    
    

    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    problem_decomposition.par_objfun.d = problem_decomposition.par_objfun.d_belief;
    problem_decomposition.par_objfun.u = problem_decomposition.par_objfun.u_belief;
    problem_decomposition.par_objfun.lb_d = problem_decomposition.lb;
    problem_decomposition.par_objfun.ub_d = problem_decomposition.ub;
    
    problem_decomposition.par_objfun.maximise = 1;    
    
    
%     problem_inner.par_objfun.problem_par_objfun{1, 1}.fix = in.par_objfun{1, 1}.fix;
    
    [ u_max_to_opt_all, fmax_to_opt_all, ~ , ~ ] = algo_inner.optimise(problem_decomposition,algo_inner.par);
    
    % take the best solution between all the population of MP AIDEA
    [fmax_to_opt, pop_max] = min(fmax_to_opt_all);
    u_max_to_opt = u_max_to_opt_all(pop_max,:);


    decomposition{1,i-in.num_functions}.FocalElement{1,ii}.upper_f = -fmax_to_opt;
    decomposition{1,i-in.num_functions}.FocalElement{1,ii}.upper_u = u_max_to_opt;
    
end


%--------------------------------------------------------------------------
% MIN
%--------------------------------------------------------------------------
if  in.flag_output.Plausibility % in.output == 1 || in.output == 2
    
    problem_decomposition.objfun = @mask_objfun_min_decomposition;
    
    
    % objective and constraints are defined in different functions
    % Function to optimise
%     fitnessfcn.obj       = problem_inner.objfun;
    % Function of constraints
    if isempty(problem_decomposition.par_objfun.constraint{1})
        problem_decomposition.constraint    = [];
    else
        problem_decomposition.constraint    = @mask_constraint_max_decomposition;
    end
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    problem_decomposition.par_objfun.d = problem_decomposition.par_objfun.d_plausibility;
    problem_decomposition.par_objfun.u = problem_decomposition.par_objfun.u_plausibility;
    problem_decomposition.par_objfun.lb_d = problem_decomposition.lb;
    problem_decomposition.par_objfun.ub_d = problem_decomposition.ub;
    
    problem_decomposition.par_objfun.maximise = 0;
    
    
    problem_decomposition.par_objfun.problem_par_objfun{1, 1}.fix = in.par_objfun{1, 1}.fix;
    
    [ u_min_to_opt_all, fmin_to_opt_all, ~ , ~ ] = algo_inner.optimise(problem_decomposition,algo_inner.par);
    

    % take the best solution between all the population of MP AIDEA
    [fmin_to_opt, pop_min] = min(fmin_to_opt_all);
    u_min_to_opt = u_min_to_opt_all(pop_min,:);

    decomposition{1,i-in.num_functions}.FocalElement{1,ii}.downer_u = u_min_to_opt;
    decomposition{1,i-in.num_functions}.FocalElement{1,ii}.downer_f = fmin_to_opt;
    
    
end




global num_maximization_decomposition
num_maximization_decomposition = num_maximization_decomposition +1;

end