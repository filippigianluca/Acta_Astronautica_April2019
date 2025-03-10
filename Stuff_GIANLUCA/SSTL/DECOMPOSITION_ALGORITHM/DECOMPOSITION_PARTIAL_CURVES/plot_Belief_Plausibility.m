function [Partial_curve] = plot_Belief_Plausibility(i, decomposition, problem_decomposition, in, problem_inner, minmax, minmin, n_obj, n_point, Partial_curve)



if in.output == 0 || in.output == 2
    
    % BELIEF
    
    for k = 1 : number_Focal_Element(problem_decomposition, problem_inner)     % vector of all the f-values of the FEs
        fmax(k) = decomposition{1,i-in.num_functions}.FocalElement{1,k}.upper_f;
    end
    
    [fmax_sort, idx_max] = sort(fmax(:));
    
    count = 0;
    bel_comulative = 0;
    for fe = idx_max'
        bel_comulative = bel_comulative+decomposition{1,i-in.num_functions}.FocalElement{1,fe}.bpa;
        count = count + 1;
        belief_partial(count) = bel_comulative;
    end
    
    
    Partial_curve{i - in.num_functions}.Belief_FE_function = fmax_sort;
    Partial_curve{i - in.num_functions}.Belief_FE_belief_partial = belief_partial;
    Partial_curve{i - in.num_functions}.Belief_FE_number = idx_max;
    
    
end


if in.output == 1 || in.output == 2
    
    % PLAUSIBILITY
    
    
    for k = 1 : number_Focal_Element(problem_decomposition, problem_inner)     % vector of all the f-values of the FEs
        
        fmin(k) = decomposition{1,i-in.num_functions}.FocalElement{1,k}.downer_f;
    end
    
    [fmin_sort, idx_min] = sort(fmin(:));
    
    count = 0;
    plausibility_comulative = 0;
    for fe = idx_min'
        plausibility_comulative = plausibility_comulative+decomposition{1,i-in.num_functions}.FocalElement{1,fe}.bpa;
        count = count + 1;
        plausibility_partial(count) = plausibility_comulative;
    end
    
    
    
    Partial_curve{i - in.num_functions}.Plausibility_FE_function = fmin_sort;
    Partial_curve{i - in.num_functions}.Plausibility_FE_plausibility_partial = plausibility_partial;
    Partial_curve{i - in.num_functions}.Plausibility_FE_number = idx_min;
    
    
end

if  in.output ~= 0 && in.output ~= 1 && in.output ~= 2
    printf('error')
    
end



end

