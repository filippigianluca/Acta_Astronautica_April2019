function [decomposition_end] = evaluate_max_min_FE_decoupled(num_function, num_FE, in, position_FE, decomposition_end, algo_inner, problem_inner, num_sample, u_max_tot, u_max_tot_Plausibility)


bpa = 1;
for k = 1:problem_inner.dim
    bpa = bpa*problem_inner.bpa{k,1}(position_FE(k));  % problem_inner.bpa{1,1}{k,1}(position_FE(k));
end

decomposition_end{num_sample}.step_one{num_function}{num_FE}.bpa = bpa;


decomposition_end{num_sample}.step_one{num_function}{num_FE}.n_FE = num_FE;
decomposition_end{num_sample}.step_one{num_function}{num_FE}.lb = problem_inner.lb;
decomposition_end{num_sample}.step_one{num_function}{num_FE}.ub = problem_inner.ub;


if in.output == 0 || in.output == 2
    
    %%% BELIEF
    problem_inner.objfun = @mask_objfun_max_decomposition;
    %
    problem_inner.par_objfun.u_belief = u_max_tot;
    
    [ u_max_to_opt, fmax_to_opt , ~ , ~ ] = algo_inner.optimise(problem_inner,algo_inner.par);
    
    
    
    decomposition_end{num_sample}.step_one{num_function}{num_FE}.upper_u = u_max_to_opt;
    decomposition_end{num_sample}.step_one{num_function}{num_FE}.upper_f = -fmax_to_opt;
    
end

if in.output == 1 || in.output == 2
    
    %%% PLAUSIBILITY
    problem_inner.objfun = @mask_objfun_min_decomposition;
    %
    problem_inner.par_objfun.u_plausibility = u_max_tot_Plausibility;
    
    [ u_min_to_opt, fmin_to_opt , ~ , ~ ] = algo_inner.optimise(problem_inner,algo_inner.par);
    
    decomposition_end{num_sample}.step_one{num_function}{num_FE}.downer_u = u_min_to_opt;
    decomposition_end{num_sample}.step_one{num_function}{num_FE}.downer_f = fmin_to_opt;
    
end

if  in.output ~= 0 && in.output ~= 1 && in.output ~= 2
    print('error')
    
end

global num_maximization_decomposition
num_maximization_decomposition = num_maximization_decomposition +1;

end