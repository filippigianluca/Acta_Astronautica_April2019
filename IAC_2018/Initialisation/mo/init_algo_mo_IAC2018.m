% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [algo_minmax, algo_outer, algo_inner] = init_algo_mo_IAC2018(minmax_problem)

%% META ALGORITHM: MACSMINMAX
algo_minmax.optimise = @optimise_mo_ausa;            % algorithm in the form [x,f,exitflag,output] = algo(problem,algo_outer,algo_inner,par_minmax)
% parameters
    if isfield(minmax_problem,'maxnfeval')
        par_minmax.maxnfeval = minmax_problem.maxnfeval;% TOTAL function evaluation limit
    else
        error('max nfeval not supplied')
    end

    par_minmax.n_d0 = 8;
    
    par_minmax.local_search_flags.validation = false;                               % Run local search (true) or function evaluation (false) in validation
    par_minmax.local_search_flags.inner = repmat(false,[1,minmax_problem.n_obj]);   % Idem inner (SO) problem
    par_minmax.local_search_flags.outer = false;                                    % Idem outer (MO) problem
    % par_minmax.tconv = 1;                               % tconv
    % par_minmax.tspr = 1;                                % tspr
algo_minmax.par_minmax = par_minmax;

if (minmax_problem.n_obj == 1)
    error('using MO settings for a SO problem!')
else
    %% ALGORITHM OUTER LOOP with archive cross-checks
    algo_outer.optimise_au = @optimise_macs_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        par_macs = struct;
        par_macs.maxnfeval = 100*minmax_problem.dim_d;      % Max function evaluations
        par_macs.popsize = 10 + minmax_problem.dim_d;       % Population size
        par_macs.rhoini = 1;                                % Initial local hypercube size
        par_macs.F = 0.9;                                   % F
        par_macs.CR = 0.9;                                  % CR
        par_macs.p_social = 1;                            % Ratio between elite and total population
        par_macs.max_arch = 100;                            % Inner archive size
        par_macs.max_arch_out = 7;                          % Output archive size
        par_macs.coord_ratio = 1;                           % Quota of coordinates examined for each set of individualistic actions
        par_macs.contr_ratio = 0.5;                         % Contraction ratio
        par_macs.draw_flag = 0;                             % Print itarations status  
        par_macs.cp = 0;                                    % constraint flag
        par_macs.MBHflag = 0;                               % number of MBH steps
        par_macs.explore_DE_strategy = 'rand';
        par_macs.social_DE_strategy = 'DE/current-to-rand/1';
        par_macs.v = 0;
        par_macs.dyn_pat_search = 0;
        par_macs.upd_subproblems = 0;
        par_macs.max_rho_contr = 5;
        par_macs.pat_search_strategy = 'standard';
        par_macs.optimal_control = 0;
        par_macs.vars_to_opt = ones(minmax_problem.dim_d,1);
        
        par_macs.use = true;
        if(par_macs.maxnfeval == 0 || par_macs.max_arch_out == 0)
            par_macs.use = false;
        end
    algo_outer.par_au = par_macs;
    

    
    %% ALGORITHM OUTER LOOP with surrogate
    algo_outer.optimise_sa = @optimise_macs_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        par_macs = struct;
        par_macs.maxnfeval = 0;%1000*minmax_problem.dim_d;      % Max function evaluations
        par_macs.popsize = 10 + minmax_problem.dim_d;       % Population size
        par_macs.rhoini = 1;                                % Initial local hypercube size
        par_macs.F = 0.9;                                   % F
        par_macs.CR = 0.9;                                  % CR
        par_macs.p_social = 0.2;                            % Ratio between elite and total population
        par_macs.max_arch = 0;                            % Inner archive size
        par_macs.max_arch_out = 0;                          % Output archive size
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
        par_macs.vars_to_opt = ones(minmax_problem.dim_d,1);
        
        par_macs.use = true;
        if(par_macs.maxnfeval == 0 || par_macs.max_arch_out == 0)
            par_macs.use = false;
        end

    algo_outer.par_sa = par_macs;

end

algo_inner = cell(1,minmax_problem.n_obj);
for obj = 1:minmax_problem.n_obj
    %% ALGORITHM INNER LOOP: MPAIDEA
    algo_inner{obj}.optimise = @optimise_mpaidea_wrapper;                % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    par_mpaidea = struct;
        par_mpaidea.nFeValMax = 1000*minmax_problem.dim_d;     % number of function evaluations for IDEA     100;%
        par_mpaidea.n_populations = 1;                         % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(5,minmax_problem.dim_d);    % number of agents in one population        
        par_mpaidea.max_LR = 10;                               % max.number of local restarts
        par_mpaidea.F = 0.75;                                    % F
        par_mpaidea.CR = 0.8;                                   % CR
        par_mpaidea.parallel = 0;                              % Parallel evaluation of individuals of the DE?        
        par_mpaidea.plots = 0;                                 % Display plots during run?
        par_mpaidea.text = 0;                                  % Display text during run?
        par_mpaidea.fmincon_set = 'interior-point';                       % choose algorithm for fmincon: 'sqp', 'interior-point', 
        %------------------------------------------------------------------
        % default parameters
        %------------------------------------------------------------------
        par_mpaidea.population=[];                              % initial population, same for all execution. Leave empty to randomize.
        par_mpaidea.DE_strategy = 1;                            % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
        par_mpaidea.prob_DE_strategy = 0.5;                     % probability or not using DE/Rand
        par_mpaidea.delta_local = 0.1;                          % dimension of the bubble for the local restart of population.
        par_mpaidea.delta_global = 0.1;                         % characteristic dimension for the global restart of the population
        par_mpaidea.rho = 0.2;                                  % contraction threshold for the population
        par_mpaidea.dd_CRF = 3;                                 % parameter for the adaptation of CRF.
        par_mpaidea.record = [1];                               % Fraction(s) of nfevalmax at which the result is saved       
        par_mpaidea.save_pop_DE = 0;                            % Save results of DE to file?
        par_mpaidea.name_save_pop_DE = 'pop_DE_';               % If yes, choose prefix of name for files:
        par_mpaidea.save_local_search = 0;                      % Save results of local search to file?
        par_mpaidea.name_save_local_search = 'minima_fmincon_'; % If yes, choose prefix for name for file:
        par_mpaidea.save_pop_LR = 0;                            % Save populations at local restart (each one saved on a different file)?
        par_mpaidea.name_save_pop_LR = 'pop_LR_';               % If yes, choose prefix of name for files:
        par_mpaidea.save_pop_GR = 0;                            % Save populations at global restart (each one saved on a different file)?
        par_mpaidea.name_save_pop_GR = 'pop_GR_';               % If yes, choose prefix of name for files
        
        par_mpaidea.fitnessfcn.obj_constr = 1;                  % Flag to 1 if objective and constraints are in the same function
        par_mpaidea.fitnessfcn.weighted   = 0;                  % How to handle constraints: set to 1 for weighted constraints with fixed weights, or to 0 for penalty with no weights
        par_mpaidea.fitnessfcn.ceq_eps    = 1e-6;               % If the constraints are handled without weights, then define a tolerance for the violation of the equality constraints
        par_mpaidea.fitnessfcn.w_ceq      = 1000;               % Weights for penalty
        par_mpaidea.fitnessfcn.w_c        = 100;
        %------------------------------------------------------------------
        
    algo_inner{obj}.par = par_mpaidea;
    

end

return