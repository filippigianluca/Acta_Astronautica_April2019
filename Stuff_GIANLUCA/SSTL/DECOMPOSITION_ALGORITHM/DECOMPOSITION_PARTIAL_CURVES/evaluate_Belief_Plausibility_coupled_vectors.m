function [decomposition, Partial_curve] = evaluate_Belief_Plausibility_coupled_vectors(in, problem_decomposition, minmax, minmin, n_obj, n_point, algo_decomposition)
% PARTIAL CURVES (EXCHANGE VARIABLES)
%
% fix d and the uncoupled vectors;
% for each coupled vector, for all the FE:
% do a maximization (minimization) over that coupled vector
%
% min_max_func_decomposizion -> evaluate max (min) and save in 'decomposition'
%                               all the informations of the FE
% plot_Belief_Plausibility   -> evaluate partial Belief curve



problem_inner = build_problem_inner(problem_decomposition, minmax, minmin, n_obj);        % fix minmax-output as decomposition-input


decomposition = cell(1,in.num_functions/2*(in.num_functions-1));
Partial_curve = cell(1,in.num_functions/2*(in.num_functions-1));


for i = in.num_functions +1 : length(in.dim_u)           % do maximization (and/or minimization) for each coupled vector
    
    
    
    if in.dim_u(i)>0                                     % DO only if there are coupled variables
        
        
        problem_inner.par_objfun.vars_to_opt = var2opt(i,problem_decomposition);      % choose the components to optimize in u
        problem_inner.dim = length(problem_inner.par_objfun.vars_to_opt);
        vars_to_opt = problem_inner.par_objfun.vars_to_opt;
        
        problem_max = problem_decomposition;
        
        problem_max.dim_u = problem_inner.dim;
        
        for index =1:problem_inner.dim
            lb_u(index,1) = {in.lb_u{n_obj}{vars_to_opt(index)}};
            ub_u(index,1) = {in.ub_u{n_obj}{vars_to_opt(index)}};
            bpa(index,1) = {in.bpa{n_obj}{vars_to_opt(index)}};
        end
        problem_max.lb_u = lb_u;
        problem_max.ub_u = ub_u;
        
        problem_max.bpa = bpa;
        
        
        
        
        
        
        
        problem_inner.bpa = problem_max.bpa;
        
        problem_decomposition.dim_u = problem_inner.dim;
        num_FE_tot = number_Focal_Element(problem_decomposition, problem_inner);
        
        for index_num_FE = 1:num_FE_tot        % for each focal element
            
            problem_inner.lb = zeros(1,problem_inner.dim);
            problem_inner.ub = zeros(1,problem_inner.dim);
            
            position_FE = position(index_num_FE, problem_inner, problem_max);
            
            for index_num_interval = 1 : problem_inner.dim    % find the domain for the chosen focal element
                
                problem_inner.lb(index_num_interval) = problem_max.lb_u{index_num_interval,1}(position_FE(index_num_interval));
                problem_inner.ub(index_num_interval) = problem_max.ub_u{index_num_interval,1}(position_FE(index_num_interval));
                
            end
            
            [decomposition] = max_min_func_decomposition(i, index_num_FE, in, position_FE, decomposition, algo_decomposition, problem_inner);
            
        end
        
        
        
        [Partial_curve] = plot_Belief_Plausibility(i, decomposition, problem_decomposition, in, problem_inner, minmax, minmin, n_obj, n_point, Partial_curve);
        
    end
    
end


% plot of the partial curves
for i = in.num_functions +1 : length(in.dim_u)
    
    if in.dim_u(i)>0
        
         
        
        if in.output == 0
            Belief_F_old = Partial_curve{i - in.num_functions}.Belief_FE_function;
            Belief_F = [Belief_F_old(1) Belief_F_old'];
            Belief = [0 Partial_curve{i - in.num_functions}.Belief_FE_belief_partial];
            
            figure
            stairs(Belief_F, Belief,'linewidth',2)
            xlabel('F')
            ylabel('partial Belief')
            title('Exchange variables (decomposition approach)')
        elseif in.output == 1
            Plausibility_F_old = Partial_curve{i - in.num_functions}.Plausibility_FE_function;
            Plausibility_F = [Plausibility_F_old(1) Plausibility_F_old'];
            Plausibility = [0 Partial_curve{i - in.num_functions}.Plausibility_FE_plausibility_partial];
            
            figure
            stairs(Plausibility_F, Plausibility,'linewidth',2)
            xlabel('F')
            ylabel('partial Plausibility')
            title('Exchange variables (decomposition approach)')
        elseif in.output == 2
            Belief_F_old = Partial_curve{i - in.num_functions}.Belief_FE_function;
            Belief_F = [Belief_F_old(1) Belief_F_old'];
            Belief = [0 Partial_curve{i - in.num_functions}.Belief_FE_belief_partial];
            Plausibility_F_old = Partial_curve{i - in.num_functions}.Plausibility_FE_function;
            Plausibility_F = [Plausibility_F_old(1) Plausibility_F_old'];
            Plausibility = [0 Partial_curve{i - in.num_functions}.Plausibility_FE_plausibility_partial];
            figure
            stairs(Belief_F, Belief,'linewidth',2)
            hold on
            stairs(Plausibility_F, Plausibility,'linewidth',2)
            xlabel('F')
            ylabel('partial Belief & Plausibility')
            title('Exchange variables (decomposition approach)')
        end
        
    end
    
end

end