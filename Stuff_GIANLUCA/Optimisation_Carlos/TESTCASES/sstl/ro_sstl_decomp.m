addpath('/home/carlos/phd/MACSMINMAX_testbench/decomposition/SSTL/Modelli_SSTL_31_01/');
addpath(genpath('/home/carlos/phd/MACSMINMAX_testbench/decomposition/SSTL/Modelli_SSTL_31_01/decomposition'));
addpath(genpath('/home/carlos/phd/MACSMINMAX_testbench/decomposition/SSTL/Modelli_SSTL_31_01/MINMIN'));
addpath(genpath('/home/carlos/phd/MACSMINMAX_testbench/decomposition/SSTL/Modelli_SSTL_31_01/sampling'));        

par_macs = struct;
        par_macs.maxnfeval = 2e4;                           % Max function evaluations
        par_macs.popsize = 10;                              % Population size
        par_macs.rhoini = 1;                                % Initial local hypercube size
        par_macs.F = 0.9;                                   % F
        par_macs.CR = 0.9;                                  % CR
        par_macs.p_social = 0.2;                            % Ratio between elite and total population
        par_macs.max_arch = 5000;                           % Inner archive size
        par_macs.max_arch_out = 100;                        % Output archive size
        par_macs.coord_ratio = 1;                           % Quota of coordinates examined for each set of individualistic actions
        par_macs.contr_ratio = 0.5;                         % Contraction ratio
        par_macs.draw_flag = 0;                             % Print itarations status  
        par_macs.cp = 0;                                    % constraint flag
        par_macs.MBHflag = 0;                               % number of MBH steps
        par_macs.explore_DE_strategy = 'rand';
        par_macs.social_DE_strategy = 'DE/current-to-rand/1';
        par_macs.v = 0;
        par_macs.dyn_pat_search = 1;
        par_macs.upd_subproblems = 0;
        par_macs.max_rho_contr = 5;
        par_macs.pat_search_strategy = 'standard';
        par_macs.optimal_control = 0;
        par_macs.vars_to_opt = [1,1];

        LB = [0.0 , 0.150];
        UB = [1.0 , 0.255];

        global iter_count_ro;
        iter_count_ro =0;
[x,fval,exitflag,output] = optimise_macs(@Fsstl3_decomp,LB,UB, par_macs);
save(strcat('ro_sstl_p3decomp_',num2str(par_macs.maxnfeval),'feval_newcost'));