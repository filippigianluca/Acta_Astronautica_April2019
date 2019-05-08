function [decomposition_TEST_FocalElement] = max_func_decomposizion_TEST(in, ii, position_FE, decomposition_TEST_FocalElement, algo_inner, problem_inner)

global num_maximization_Belief_exact;
num_maximization_Belief_exact = num_maximization_Belief_exact + 1;

%% MAX
if in.output == 0 || in.output == 2
    problem_inner.objfun = @mask_objfun_max_decomposition;
    
    [ u_max_to_opt, fmax_to_opt , ~ , ~ ] = algo_inner.optimise(problem_inner,algo_inner.par);
    
    
    
    decomposition_TEST_FocalElement{1,ii}.upper_u = u_max_to_opt;
    decomposition_TEST_FocalElement{1,ii}.upper_f = -fmax_to_opt;
    
end

%% MIN
if in.output == 1 || in.output == 2
    
    problem_inner.objfun = @mask_objfun_min_decomposition;
    
    [ u_min_to_opt, fmin_to_opt , ~ , ~ ] = algo_inner.optimise(problem_inner,algo_inner.par);
    
    decomposition_TEST_FocalElement{1,ii}.downer_u = u_min_to_opt;
    decomposition_TEST_FocalElement{1,ii}.downer_f = fmin_to_opt;
end



bpa = 1;
for k = 1:problem_inner.dim
    bpa = bpa*problem_inner.bpa{k,1}(position_FE(k));  % problem_inner.bpa{1,1}{k,1}(position_FE(k));
end

decomposition_TEST_FocalElement{1,ii}.bpa = bpa;
decomposition_TEST_FocalElement{1,ii}.n_FE = ii;
decomposition_TEST_FocalElement{1,ii}.lb = problem_inner.lb;
decomposition_TEST_FocalElement{1,ii}.ub = problem_inner.ub;


end