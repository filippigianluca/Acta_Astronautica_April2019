% test monotonicity
clear all; close all; clc



load('partial_curve_10000feval.mat')

N_val = 30000;
%--------------------------------------------------------------------------
options.u_fix = minmax.u;

% Function to optimise ----------------------------------------------------
fitnessfcn.obj = @mask_test_monotonicity;
fitnessfcn.constr = [] ;
fitnessfcn.obj_constr = 0;
fitnessfcn.weighted = 0;
%--------------------------------------------------------------------------

LB_tot = problem.lb_u{1, 1};
UB_tot = problem.ub_u{1, 1};  
for k = 1:length(Partial_curve)    
    if ~isempty(Partial_curve{k})
        
        u_fixed_Partial = problem.order_dim_u(sum(problem.dim_u_i(1:problem.num_functions+k-1))+1:sum(problem.dim_u_i(1:problem.num_functions+k))); 
        
        u_opti = [];
        for i = 1:length(problem.order_dim_u)
            if ~any(u_fixed_Partial == i)
                u_opti = [u_opti i];
            end
        end
%         u_opti = [1:u_fixed_Partial(1)-1  u_fixed_Partial(end)+1:length(problem.order_dim_u)];
        
        problem.lb = [];
        problem.ub = [];
        for m=u_opti
           problem.lb = [problem.lb  LB_tot{m}(1)];
           problem.ub = [problem.ub  UB_tot{m}(2)];
        end
        
        problem.dim = length(problem.lb);
        problem.fitnessfcn = fitnessfcn;
        problem.objfun = @IAC_objfun;
        
        problem.par_objfun{1, 1}.u_fix = minmax.u;
        problem.par_objfun{1, 1}.d_fix = minmax.d;
        problem.par_objfun{1, 1}.u_fixed_Partial = u_fixed_Partial;
        problem.par_objfun{1, 1}.u_opti = u_opti;
        
        par = algo_decomposition.par;
        par.nFeValMax = N_val;
        
        
        algo_decomposition.par.fmincon_set = 'sqp';
        for i=1:length(Partial_curve{1, k}.u_coupled)
            
          problem.par_objfun{1, 1}.u_coupled = Partial_curve{1, k}.u_coupled(i, :); 
          
          [ x, fval, exitflag, output ] = optimise_mpaidea_wrapper(problem,par);
            
          [fval_max, pos_fval_max] = max(-fval);
          
          F_MAX(k, i) = fval_max;
          X_MAX{k, i} = x(pos_fval_max, :);
        end
    end
end

save(strcat('solution_max_ui_uij_fixed_',num2str(N_val),'feval'))


