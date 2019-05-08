% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [decomposition, Partial_curve] = evaluate_Belief_Plausibility_coupled_vectors(problem, minmax, minmin, n_obj, n_point, algo)

% PARTIAL CURVES (EXCHANGE VARIABLES)
%
% fix d and the uncoupled vectors;
% for each coupled vector, for all the FE:
% do a maximization (minimization) over that coupled vector
%
% min_max_func_decomposizion -> evaluate max (min) and save in 'decomposition'
%                               all the informations of the FE
% plot_Belief_Plausibility   -> evaluate partial Belief curve

disp(strcat('Running Decomposition Coupled Belief...'));



%% initialisation
decomposition = cell(1,problem.num_functions/2*(problem.num_functions-1));
Partial_curve = cell(1,problem.num_functions/2*(problem.num_functions-1));

problem_coupled = build_problem_inner(problem, minmax, minmin, n_obj);        % fix minmax-output as decomposition-input




%% evaluate Belief and/or Plausibility curve(s) for each coupled vector
for i = problem.num_functions +1 : length(problem.dim_u_i)          
    
    clearvars lb_u  ub_u bpa pos_FE position_FE
    
    vars_to_opt = var2opt(i,problem);      % choose the components to optimize in u in the input order
    
    % DO only if there are coupled variables
    if problem.dim_u_i(i)>0                                     &&  ...
            ( problem.num_samples{i-problem.num_functions}~=1   ||  ...
             length(problem.sample.FE_uSample_M{problem.indexing.matrix2linear{i}(1), problem.indexing.matrix2linear{i}(2)})>1)

           
        
        problem_coupled.dim = length(vars_to_opt);
        
        problem_coupled.par_objfun.vars_to_opt = vars_to_opt;
        
        
        problem_max = problem;
        
        problem_max.dim_u = problem_coupled.dim;
        
        for index =1:problem_coupled.dim
            lb_u(index,1) = {problem.lb_u{n_obj}{vars_to_opt(index)}};
            ub_u(index,1) = {problem.ub_u{n_obj}{vars_to_opt(index)}};
            bpa(index,1) = {problem.bpa{n_obj}{vars_to_opt(index)}};
            %             lb_u(index,1) =
            %             {in.lb_u{n_obj}{problem_decomposition.order_dim_u(vars_to_opt(index))}}; 
            %             ub_u(index,1) = {in.ub_u{n_obj}{problem_decomposition.order_dim_u(vars_to_opt(index))}};
            %             bpa(index,1)  = {in.bpa{n_obj}{problem_decomposition.order_dim_u(vars_to_opt(index))}};
            lengt_bpa(1, index) = length(bpa{index});
        end
        problem_max.lb_u = lb_u;
        problem_max.ub_u = ub_u;
        
        problem_max.bpa = bpa;
        
        
        
        %------------------------------------------------------------------
        % initialise "pos_FE" to find the position in FE space with "ind2sub" 
        %------------------------------------------------------------------    
        pos_FE = '[position_FE(1)';
        for iii=2:length(lb_u)
            pos_FE = strcat(pos_FE, ',position_FE(',num2str(iii),')');
        end
        pos_FE = strcat(pos_FE, ']');
        %------------------------------------------------------------------

        
        problem_coupled.bpa = problem_max.bpa;
        
        problem.dim_u = problem_coupled.dim;
        num_FE_tot = number_Focal_Element(problem, problem_coupled);
        
        for index_num_FE = 1:num_FE_tot        % for each focal element
            
            problem_coupled.lb = zeros(1,problem_coupled.dim);
            problem_coupled.ub = zeros(1,problem_coupled.dim);
            
            
            %             position_FE = position(index_num_FE, problem_inner, problem_max);
            eval(strcat(pos_FE, ' = ind2sub([',num2str(lengt_bpa),'],', num2str(index_num_FE), ');'));
            
            for index_num_interval = 1 : problem_coupled.dim    % find the domain for the chosen focal element
                
                problem_coupled.lb(index_num_interval) = problem_max.lb_u{index_num_interval,1}(position_FE(index_num_interval));
                problem_coupled.ub(index_num_interval) = problem_max.ub_u{index_num_interval,1}(position_FE(index_num_interval));
                
            end
            
            [decomposition] = max_min_func_decomposition(i, index_num_FE, problem, position_FE, decomposition, algo, problem_coupled);
            
            
        end
        
        
        
        [Partial_curve] = plot_Belief_Plausibility(i, decomposition, problem, problem, problem_coupled, minmax, minmin, n_obj, n_point, Partial_curve);
        
    end
    
end




end