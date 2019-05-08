function [decomposition_end] = evaluate_combination_FE_decoupled(decomposition_end, num_sample, problem_decomposition, minmax, minmin, u_max_tot, Sample, PS, u_max_tot_Plausibility, n_obj, n_point, in)


%% FIND number of combinations for the final curve
num_combinations_max = 1;
num_functions_tot = length(decomposition_end{num_sample}.step_one);

num_functions_tot_new = 0;
for num_functions = 1:num_functions_tot
    if isempty(decomposition_end{num_sample}.step_one{num_functions})==0
        num_combinations_max = num_combinations_max*length(decomposition_end{num_sample}.step_one{num_functions});
        
        num_functions_tot_new = num_functions_tot_new + 1;
        sub_not_empty(num_functions_tot_new) = num_functions;
    end
end
%num_functions_tot = num_functions_tot_new;

%% BELIEF

if in.output == 0 || in.output == 2
    
    d_belief = minmax.d;
    f_cost = problem_decomposition.objfun{n_obj}(d_belief, u_max_tot, problem_decomposition);
    
    
    
% EVALUATE bpa-scale: for each sample a scaled curve is evaluated. 

    bpa_scale = 1;
    
    for j = 1:length(PS)
        if in.dim_u(in.num_functions + j)>0
            if PS(j) == 1
                bpa_scale = bpa_scale*Sample{j}.Belief_f(PS(j));
                
            else
                bpa_scale = bpa_scale*(Sample{j}.Belief_f(PS(j)) - Sample{j}.Belief_f(PS(j)-1));
                
            end
        end
    end
      
    
% RECONSTRUCTION of all the focal elements of the final scaled curves 

    for num_combination = 1:num_combinations_max
        
        [P_FE] = position_sum_function(num_combinations_max, num_functions_tot, num_combination, decomposition_end, num_sample);
        
        bpa =1;
        F = 0;
        
        for i = 1:length(P_FE)
            bpa = bpa*decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.bpa;
            F = F + decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.upper_f;
            
        end
        
        
        decomposition_end{num_sample}.step_two{num_combination}.bpa = bpa*bpa_scale;
        
        N=0;
        for i=1:num_functions_tot
            if in.dim_u(i)>0
                N = N+1;
            end
        end
        decomposition_end{num_sample}.step_two{num_combination}.F = F - (N-1)*f_cost;
       
    end
    
    
end

%% PLAUSIBILITY

if in.output == 1 || in.output == 2
    
    d_plausibility = minmin.d;
    f_cost_Plausibility = problem_decomposition.objfun{n_obj}(d_plausibility, u_max_tot_Plausibility, problem_decomposition);
    
    
% EVALUATE bpa-scale: for each sample a scaled curve is evaluated. 
    
    bpa_scale_Plausibility = 1;
    for j = 1:length(PS)
        if in.dim_u(in.num_functions + j)>0
            if PS(j) == 1
                
                bpa_scale_Plausibility = bpa_scale_Plausibility*Sample{j}.Plausibility_f(PS(j));
            else
                
                bpa_scale_Plausibility = bpa_scale_Plausibility*(Sample{j}.Plausibility_f(PS(j)) - Sample{j}.Plausibility_f(PS(j)-1));
            end
        end
    end
    
    
% RECONSTRUCTION of all the focal elements of the final scaled curves   

    for num_combination = 1:num_combinations_max
        
        [P_FE] = position_sum_function(num_combinations_max, num_functions_tot, num_combination, decomposition_end, num_sample);
        
        bpa =1;
        
        F_Plausibility = 0;
        for i = 1:length(P_FE)
            bpa = bpa*decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.bpa;
            
            F_Plausibility = F_Plausibility + decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.downer_f;
        end
        
        
        
        decomposition_end{num_sample}.step_two{num_combination}.bpa_Plausibility = bpa*bpa_scale_Plausibility;
        
        
        
        decomposition_end{num_sample}.step_two{num_combination}.F_Plausibility = F_Plausibility - (num_functions_tot-1)*f_cost_Plausibility;
    end
    
    
end



end