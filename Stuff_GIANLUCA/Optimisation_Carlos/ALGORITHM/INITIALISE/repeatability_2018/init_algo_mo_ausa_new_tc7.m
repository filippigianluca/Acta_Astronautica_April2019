function [algo_minmax, algo_outer, algo_inner] = init_algo_mo_ausa_new_tc7(minmax_problem);

%% META ALGORITHM: MACSMINMAX
algo_minmax.optimise = @optimise_ausa_trisurr;            % algorithm in the form [x,f,exitflag,output] = algo(problem,algo_outer,algo_inner,par_minmax)
% parameters
    if isfield(minmax_problem,'maxnfeval')
        par_minmax.maxnfeval = 7e5;% TOTAL function evaluation limit
    else
        error('max nfeval not supplied')
    end
    par_minmax.max_iter = inf;

    par_minmax.verbosity = true;
    % initial sample of D. a maximisation will be run for each
    par_minmax.n_d0 = 8;

    % Cross-check options
    cc_opt = struct;
        cc_opt.to_crosscheck = 1;             % 0 cross-checks no one, 1 cross-checks the reference, 2 cross-checks everyone in the archive
        % cc_opt.max_dom_idx_ref = 5;           % guys with dom_idx <= this are considered reference. Recommended >> 0 for SO and lower for MO. [NOT IMPLEMENTED YET]
        % cc_opt.max_ref_size = inf;            [NOT IMPLEMENTED]
        cc_opt.ls_flag_crosscheck = 1;        % 0 for no local search, 1 for local search from the argmax. in Au, 2 for local search at every cross-check
        cc_opt.restrain_u_record = true;      % if true, will purge from u_record the guys that no longer maximise for a d in d_record
        cc_opt.tol = 1e8;                     % tolerance for considering two (scaled) d's or u's equal. gets overriden by choice in optimise_blah                
        cc_opt.dim_d = minmax_problem.dim_d;
        cc_opt.dim_u = minmax_problem.dim_u;
        cc_opt.n_obj = minmax_problem.n_obj;
        cc_opt.algo_ls = cell(1,cc_opt.n_obj);
            for obj = 1:cc_opt.n_obj
                %% ALGORITHM INNER LOOP: FMINCON's SQP
                cc_opt.algo_ls{obj}.optimise = @optimise_fmincon_wrapper;                % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
                par_cc_ls = struct;
                
                % default parameters --> random choice of u0 // the u0 will be chosen by the cross-checks
                cc_opt.algo_ls{obj}.par = par_cc_ls;
            end

    par_minmax.cc_opt = cc_opt;

    par_minmax.alternate_outer = false;   % if alternate, will alternate in outer loop the methods
                                          % in order: Phi / s(d,u) subproblem / min_D (max_Au) (if .use)
                                          % switching one to another only if the reference gets stuck.
algo_minmax.par_minmax = par_minmax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% OUTER LOOP %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S(d,u) as a min-max subproblem
algo_outer.par_sp.use = false;
algo_outer.par_sp.init_subproblem = @init_subproblem_sdu;
algo_outer.par_sp.init_subalgo = @init_subalgo_sdu_ausa;
algo_outer.par_sp.n_du0 = 10;               % initial sample of DU to maintain for S(d,u). a function evaluation for each (if par_sp.use)
algo_outer.par_sp.n_u_ref_diversity = 15;   % each iteration of the minmax loop, this number of fevals is dedicated to adding u-diversity around the reference for S(d,u) (if par_sp.use)
algo_outer.par_sp.size_d_box_ref_diversity = 1e-5;  % defines the entourage of a reference d.

% AU an Phi
if (minmax_problem.n_obj == 1)
    warning('MO_noncon_testbench?')
    %% ALGORITHM OUTER LOOP for cross-checks: MPAIDEA
    % algo_outer.optimise_au = @optimise_mpaidea_wrapper;       % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %     par_mpaidea.nFeValMax = 4.5e2*minmax_problem.dim_d;     % number of function evaluations for IDEA
    %     par_mpaidea.n_populations = 1;                        % number of populations, if no adaptive behaviour should set to 1
    %     par_mpaidea.n_agents = max(4,minmax_problem.dim_d);   % number of agents in one population
    %     par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
    %     par_mpaidea.max_LR = 5;                              % max.number of local restarts
    %     par_mpaidea.DE_strategy = 1;                          % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
    %     par_mpaidea.prob_DE_strategy = 0.5;                   % probability or not using DE/Rand
    %     par_mpaidea.delta_local = 0.1;                        % dimension of the bubble for the local restart of population.
    %     par_mpaidea.delta_global = 0.1;                       % characteristic dimension for the global restart of the population
    %     par_mpaidea.rho = 0.25;                               % contraction threshold for the population
    %     par_mpaidea.dd_CRF = 3;                               % parameter for the adaptation of CRF.
    %     par_mpaidea.F = [];                                  % F
    %     par_mpaidea.CR = [];                                 % CR
    %     par_mpaidea.text = 0;                                 % verbosity
    %     par_mpaidea.record = [1];                             % Fraction(s) of nfevalmax at which the result is saved
        
    %     par_mpaidea.use = true;
    %     if(par_mpaidea.nFeValMax <= 0)
    %         par_mpaidea.use = false;
    %     end
    % algo_outer.par_au = par_mpaidea;
    
    % algo_outer.par_au.local_search_flag = false;

    % %% ALGORITHM OUTER LOOP for SURR: MPAIDEA
    % algo_outer.optimise_sa = @optimise_mpaidea_wrapper;       % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %     par_mpaidea.nFeValMax = 2e4*minmax_problem.dim_d;     % number of function evaluations for IDEA
    %     par_mpaidea.n_populations = 1;                        % number of populations, if no adaptive behaviour should set to 1
    %     par_mpaidea.n_agents = max(4,minmax_problem.dim_d);   % number of agents in one population
    %     par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
    %     par_mpaidea.max_LR = 5;                              % max.number of local restarts
    %     par_mpaidea.DE_strategy = 1;                          % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
    %     par_mpaidea.prob_DE_strategy = 0.5;                   % probability or not using DE/Rand
    %     par_mpaidea.delta_local = 0.1;                        % dimension of the bubble for the local restart of population.
    %     par_mpaidea.delta_global = 0.1;                       % characteristic dimension for the global restart of the population
    %     par_mpaidea.rho = 0.25;                               % contraction threshold for the population
    %     par_mpaidea.dd_CRF = 3;                               % parameter for the adaptation of CRF.
    %     par_mpaidea.F = [];                                  % F
    %     par_mpaidea.CR = [];                                 % CR
    %     par_mpaidea.text = 0;                                 % verbosity
    %     par_mpaidea.record = [1];                             % Fraction(s) of nfevalmax at which the result is saved
        
    %     par_mpaidea.use = true;
    %     if(par_mpaidea.nFeValMax <= 0)
    %         par_mpaidea.use = false;
    %     end
    % algo_outer.par_sa = par_mpaidea;
else

%% ALGORITHM OUTER LOOP with archive cross-checks
    algo_outer.optimise_au = @optimise_macs_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        par_macs = struct;
        par_macs.maxnfeval = 150*minmax_problem.dim_d;      % Max function evaluations
        par_macs.popsize = 8 + minmax_problem.dim_d;       % Population size
        par_macs.rhoini = 1;                                % Initial local hypercube size
        par_macs.F = 0.9;                                   % F
        par_macs.CR = 0.9;                                  % CR
        par_macs.p_social = 0.2;                            % Ratio between elite and total population
        par_macs.max_arch = 500;                            % Inner archive size
        par_macs.max_arch_out = 5;                          % Output archive size
        par_macs.coord_ratio = 1;                           % Quota of coordinates examined for each set of individualistic actions
        par_macs.contr_ratio = 0.5;                         % Contraction ratio
        par_macs.draw_flag = -1;                             % Print itarations status  
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
        par_macs.vars_to_opt = ones(minmax_problem.dim_d,1);
        par_macs.use = true;
        if(par_macs.maxnfeval <= 0 || par_macs.max_arch_out <= 0)
            par_macs.use = false;
        end
    algo_outer.par_au = par_macs;

    algo_outer.par_au.local_search_flag = false;

    %% ALGORITHM OUTER LOOP with surrogate
    algo_outer.optimise_sa = @optimise_macs_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        par_macs = struct;
        par_macs.maxnfeval = 1e3*minmax_problem.dim_d;      % Max function evaluations
        par_macs.popsize = 8 + minmax_problem.dim_d;        % Population size
        par_macs.rhoini = 1;                                % Initial local hypercube size
        par_macs.F = 0.9;                                   % F
        par_macs.CR = 0.9;                                  % CR
        par_macs.p_social = 0.2;                            % Ratio between elite and total population
        par_macs.max_arch = 500;                            % Inner archive size
        par_macs.max_arch_out = 5;                          % Output archive size
        par_macs.coord_ratio = 1;                           % Quota of coordinates examined for each set of individualistic actions
        par_macs.contr_ratio = 0.5;                         % Contraction ratio
        par_macs.draw_flag = -1;                             % Print itarations status  
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
        par_macs.vars_to_opt = ones(minmax_problem.dim_d,1);
        par_macs.use = true;
        if(par_macs.maxnfeval <= 0 || par_macs.max_arch_out <= 0)
            par_macs.use = false;
        end
    algo_outer.par_sa = par_macs;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% INNER LOOP %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

algo_inner = cell(1,minmax_problem.n_obj);
obj = 1;

    % %% ALGORITHM INNER LOOP: FMINCON's SQP
    % algo_inner{obj}.optimise = @optimise_fmincon_wrapper;                % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    % par_local = struct;
    %     % default parameters --> random choice of u0
    % algo_inner{obj}.par = par_local;
    % algo_inner{obj}.par.psi_surr = true;

    %% ALGORITHM INNER LOOP: MPAIDEA
    algo_inner{obj}.optimise = @optimise_mpaidea_wrapper;                % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    par_mpaidea = struct;
        par_mpaidea.nFeValMax = 2e2*minmax_problem.dim_u;     % number of function evaluations for IDEA
        par_mpaidea.n_populations = 1;                        % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(4,minmax_problem.dim_u);   % number of agents in one population
        par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
        par_mpaidea.max_LR = 1;                               % max.number of local restarts
        par_mpaidea.DE_strategy = 1;                          % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
        par_mpaidea.prob_DE_strategy = 0.5;                   % probability or not using DE/Rand
        par_mpaidea.delta_local = 0.1;                        % dimension of the bubble for the local restart of population.
        par_mpaidea.delta_global = 0.1;                       % characteristic dimension for the global restart of the population
        par_mpaidea.rho = 0.25;                               % contraction threshold for the population
        par_mpaidea.dd_CRF = 3;                               % parameter for the adaptation of CRF.
        par_mpaidea.F = 0.9;%[];                                  % F
        par_mpaidea.CR = 0.2;%[];                                 % CR
        par_mpaidea.text = 0;    
        par_mpaidea.record = [1];                             % Fraction(s) of nfevalmax at which the result is saved
    algo_inner{obj}.par = par_mpaidea;
    algo_inner{obj}.par.psi_surr = true;

obj = 2;

    % %% ALGORITHM INNER LOOP: FMINCON's SQP
    % algo_inner{obj}.optimise = @optimise_fmincon_wrapper;                % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    % par_local = struct;
    %     % default parameters --> random choice of u0
    % algo_inner{obj}.par = par_local;
    % algo_inner{obj}.par.psi_surr = true;

    %% ALGORITHM INNER LOOP: MPAIDEA
    algo_inner{obj}.optimise = @optimise_mpaidea_wrapper;                % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    par_mpaidea = struct;
        par_mpaidea.nFeValMax = 1e3*minmax_problem.dim_u;     % number of function evaluations for IDEA
        par_mpaidea.n_populations = 1;                        % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(4,minmax_problem.dim_u);   % number of agents in one population
        par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
        par_mpaidea.max_LR = 1;                               % max.number of local restarts
        par_mpaidea.DE_strategy = 1;                          % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
        par_mpaidea.prob_DE_strategy = 0.5;                   % probability or not using DE/Rand
        par_mpaidea.delta_local = 0.1;                        % dimension of the bubble for the local restart of population.
        par_mpaidea.delta_global = 0.1;                       % characteristic dimension for the global restart of the population
        par_mpaidea.rho = 0.25;                               % contraction threshold for the population
        par_mpaidea.dd_CRF = 3;                               % parameter for the adaptation of CRF.
        par_mpaidea.F = 0.9;%[];                                  % F
        par_mpaidea.CR = 0.2;%[];                                 % CR
        par_mpaidea.text = 0;    
        par_mpaidea.record = [1];                             % Fraction(s) of nfevalmax at which the result is saved
    algo_inner{obj}.par = par_mpaidea;
    algo_inner{obj}.par.psi_surr = true;

obj = 3;

    % %% ALGORITHM INNER LOOP: FMINCON's SQP
    % algo_inner{obj}.optimise = @optimise_fmincon_wrapper;                % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    % par_local = struct;
    %     % default parameters --> random choice of u0
    % algo_inner{obj}.par = par_local;
    % algo_inner{obj}.par.psi_surr = true;

    %% ALGORITHM INNER LOOP: MPAIDEA
    algo_inner{obj}.optimise = @optimise_mpaidea_wrapper;                % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    par_mpaidea = struct;
        par_mpaidea.nFeValMax = 8e2*minmax_problem.dim_u;     % number of function evaluations for IDEA
        par_mpaidea.n_populations = 1;                        % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(4,minmax_problem.dim_u);   % number of agents in one population
        par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
        par_mpaidea.max_LR = 1;                               % max.number of local restarts
        par_mpaidea.DE_strategy = 1;                          % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
        par_mpaidea.prob_DE_strategy = 0.5;                   % probability or not using DE/Rand
        par_mpaidea.delta_local = 0.1;                        % dimension of the bubble for the local restart of population.
        par_mpaidea.delta_global = 0.1;                       % characteristic dimension for the global restart of the population
        par_mpaidea.rho = 0.25;                               % contraction threshold for the population
        par_mpaidea.dd_CRF = 3;                               % parameter for the adaptation of CRF.
        par_mpaidea.F = 0.9;%[];                                  % F
        par_mpaidea.CR = 0.2;%[];                                 % CR
        par_mpaidea.text = 0;    
        par_mpaidea.record = [1];                             % Fraction(s) of nfevalmax at which the result is saved
    algo_inner{obj}.par = par_mpaidea;
    algo_inner{obj}.par.psi_surr = true;
return