function [bel] = sstl2p_decomp(design)
bel = 0.0;
objective = 1;
% n_iter = 7;
% presample = 200;
% opt_count = 0;

init0 = str2func(strcat('init_algo_so_sstl_corners_decomp'));
% savefolder = strcat('RESULTS/SC/');

% global nfevalglobal;
% nfevalglobal = 0;

% global f_history;
% f_history=[];

%% initialise problem
[ problem_0 ] = init_problem_sstl2_decomp();
% problem_0.split_coordinates = [];

% dim_u = problem_0.dim_u;
% lb_u = problem_0.lb_u{objective};
% ub_u = problem_0.ub_u{objective};
% problem_list_next = [problem_0];
% n_int_orig = cellfun('size',problem_0.lb_u{objective},2);

[ ~, ~, algo_inner ] = init0(problem_0);
problem_max_u = build_metaproblem_macsminmax_inner(problem_0);
problem_max_u.par_objfun.objective = objective;
problem_max_u.par_objfun.d = (design-problem_0.lb_d')./(problem_0.ub_d'-problem_0.lb_d');

%% optimise
[ umax, fmax , ~ , output_aux] = algo_inner.optimise(problem_max_u,algo_inner.par);
% opt_count = opt_count+1;
fmax = -fmax;
problem_0.provisional_fmax = fmax;

if(fmax<=0)
    bel = 1.0;
else
    [ umin, fmin , ~ , output_aux] = algo_inner.minimise(problem_max_u,algo_inner.par);
    fmin = -fmin;
    
    if (fmin<0) %otherwise nothing to do but exit and return bel = 0
        % clear all; close all; clc;

% DECOMPOSITION
% addpath(genpath('/home/carlos/phd/MACSMINMAX_testbench/SSTL/Modelli_sstl_31_01/MINMAX_NOsurrogate_1obj'));

% addpath(genpath('/home/carlos/phd/MACSMINMAX_testbench/SSTL/Modelli_sstl_31_01/test'));
% addpath(genpath('/home/carlos/phd/MACSMINMAX_testbench/SSTL/Modelli_sstl_31_01/function'));
% addpath(genpath('/home/carlos/phd/MACSMINMAX_testbench/SSTL/Modelli_sstl_31_01/smart-sys-master'));
% load('D:\smart-o2c\TC-DECOMPOSITION-CARLOS\TC2_decomposition.mat')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_decomposition:                                                    
%                                                                       
% 1) function  'evaluate_minmax'                              
%                                                                         
% 2) function  'evaluate_minmin'                          
%
% 3) function  'evaluate_Belief_coupled_vectors'  
%
% 4) function  'Sampling_Belief_Plausibility'
%
% 5) function  'reconstruction_Belief_Plausibility'           
%
% 6) function  'evaluate_Belief_Plausibility_exact'
%                                                                          
% 7) plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% INITIALISATION

init = str2func('init_algo_so_sstl_corners_decomp');

[in, problem] = init_decomposition_sstl_decomp();
n_obj_tot = 1;

for n_point = 1%:n_points_tot                  % all the design points 
    
    for n_obj = 1:n_obj_tot
        
            map_info = problem_max_u.par_objfun.map_u_info{objective};
            minmax.d = design;%(design-problem_0.lb_d')./(problem_0.ub_d'-problem_0.lb_d');
            minmax.u = map_affine(umax,map_info);
            minmax.f = problem.objfun{n_obj}(minmax.d,minmax.u,problem.par_objfun{n_obj});
            
            minmin.d = design;%(design-problem_0.lb_d')./(problem_0.ub_d'-problem_0.lb_d');
            minmin.u = map_affine(umin,map_info);
            minmin.f = problem.objfun{n_obj}(minmin.d,minmin.u,problem.par_objfun{n_obj});
            %             minmax.d = ds(n_point, :);                                     % for Carlos' testcase
%             minmax.u = us_max{n_obj}(n_point,:); 
%             minmax.f = fs_max(n_point, n_obj); 
% 
%             minmin.d = ds(n_point, :); 
%             minmin.u = us_min{n_obj}(n_point,:); 
%             minmin.f = fs_min(n_point, n_obj);           
        
        %%  
        
        global num_maximization_decomposition;
        num_maximization_decomposition = 0;
        
        [decomposition, Partial_curve] = evaluate_Belief_Plausibility_coupled_vectors(in, problem, minmax, minmin, n_obj, n_point);
        
        
        %% sampling: 
        
        [Sample] = Sampling_Belief_Plausibility(Partial_curve, decomposition, in, minmax, minmin, n_obj, n_point);
        
        
        %% reconstruction
        
        [decomposition_end, Plot_scaled, Plot_decomposition, num_sample_tot] = reconstruction_Belief_Plausibility(in, problem, minmax, minmin, Sample, n_obj, n_point);
        
       
        %% validation of Belief curve: run optimization for each focal element
        
%         global num_maximization_Belief_exact;
%         num_maximization_Belief_exact = 0;
        
        
%         [decomposition_TEST, Plot_TEST] = evaluate_Belief_Plausibility_exact(in, problem, minmax, minmin, n_obj, n_point);
%         
% figure    
% hold on
% plot(Plot_TEST.f{1}, Plot_TEST.Belief{1}, 'k-.', 'linewidth', 2)
% plot(Plot_scaled{1}.f, Plot_decomposition.bl{1}, 'r-.', 'linewidth', 1)
% legend('Belief exact','Decomposition')
% hold off
        
        
        %% chart_(true and approxiated)
        
%         [Plot_end, Plot_decomposition] = plot_all(Plot_end, decomposition_end, minmax, minmin, Plot_TEST, n_obj);
        
        %create results directory
        
        % savefolder = strcat('RESULTS/SO_function_sstl/');
        % mkdir(savefolder);
        % save(strcat(savefolder,'belief'));
        
        
    end
    
end

    global num_maximization_decomposition;
        [~,bin] = histc(0.0,Plot_decomposition.f{1});
        global iter_count_ro;
        iter_count_ro = iter_count_ro +1
        bel = Plot_decomposition.bl{1}(bin)
    end
end
return


% %create results directory
% mkdir(savefolder);
% save(strcat(savefolder,'SC_',num2str(runid)));

% end