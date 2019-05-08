function [decomposition] = max_min_func_decomposition(i, ii, in, position_FE, decomposition, algo_inner, problem_inner)


bpa = 1;
for k = 1:problem_inner.dim
    bpa = bpa*problem_inner.bpa{k,1}(position_FE(k));  % product of the bpa of the coupled components 
end

decomposition{1,i-in.num_functions}.FocalElement{1,ii}.bpa = bpa;
decomposition{1,i-in.num_functions}.FocalElement{1,ii}.n_FE = ii;
decomposition{1,i-in.num_functions}.FocalElement{1,ii}.lb = problem_inner.lb;
decomposition{1,i-in.num_functions}.FocalElement{1,ii}.ub = problem_inner.ub;


if in.output == 0 || in.output == 2
    
        
    %%% MAX
    problem_inner.objfun = @mask_objfun_max_decomposition;
    
    [ u_max_to_opt, fmax_to_opt , ~ , ~ ] = algo_inner.optimise(problem_inner,algo_inner.par);
    
    decomposition{1,i-in.num_functions}.FocalElement{1,ii}.upper_u = u_max_to_opt;
    decomposition{1,i-in.num_functions}.FocalElement{1,ii}.upper_f = -fmax_to_opt;
    
end  

if in.output == 1 || in.output == 2
    
    
    %%% MIN
    problem_inner.objfun = @mask_objfun_min_decomposition;
    
    [ u_min_to_opt, fmin_to_opt , ~ , ~ ] = algo_inner.optimise(problem_inner,algo_inner.par);
    
    decomposition{1,i-in.num_functions}.FocalElement{1,ii}.downer_u = u_min_to_opt;
    decomposition{1,i-in.num_functions}.FocalElement{1,ii}.downer_f = fmin_to_opt;
 
    
end  

if  in.output ~= 0 && in.output ~= 1 && in.output ~= 2
    print('error')
    
end

global num_maximization_decomposition
    num_maximization_decomposition = num_maximization_decomposition +1;

end