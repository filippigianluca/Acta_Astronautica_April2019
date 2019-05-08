% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [algo_minmax, algo_outer, algo_inner, algo_decomposition] = init_algo_so_constr_ebro_IAC_2018_old(problem_ebro)

%% META ALGORITHM: MACSMINMAX
algo_minmax.optimise = @optimise_so;            % algorithm in the form [x,f,exitflag,output] = algo(problem,algo_outer,algo_inner,par_minmax)
% parameters
    if isfield(problem_ebro,'maxnfeval')
        par_minmax.maxnfeval = problem_ebro.maxnfeval;% TOTAL function evaluation limit
    else
        par_minmax.maxnfeval = 5e6;                     %error('max nfeval not supplied') % 100;%
    end
    par_minmax.n_d0 = 1;
    
    par_minmax.local_search_flags.validation = false;    % Run local search (true) or function evaluation (false) in validation
    par_minmax.local_search_flags.inner = false;        % Idem inner (SO) problem
    par_minmax.local_search_flags.outer = false;       % Idem outer (MO) problem

algo_minmax.par_minmax = par_minmax;


%% ALGORITHM OUTER LOOP
if (problem_ebro.n_obj == 1)
    % ALGORITHM OUTER LOOP: MPAIDEA
    algo_outer.optimise = @optimise_mpaidea_wrapper;          % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)

        par_mpaidea.nFeValMax = 15000;                          % number of function evaluations for IDEA     100;%
        par_mpaidea.n_populations = 1;                         % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(5,problem_ebro.dim_d);      % number of agents in one population        
        par_mpaidea.max_LR = 20;                               % max.number of local restarts
        par_mpaidea.F = 0.75;                                    % F
        par_mpaidea.CR = 0.8;                                   % CR
        par_mpaidea.parallel = 0;                              % Parallel evaluation of individuals of the DE?        
        par_mpaidea.plots = 0;                                 % Display plots during run?
        par_mpaidea.text = 0;                                  % Display text during run?
        par_mpaidea.fmincon_set = 'sqp';                       % choose algorithm for fmincon: 'sqp', 'interior-point', 
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
        
  algo_outer.par = par_mpaidea;


else
    warning('using SO settings for a MO problem!')
    %% ALGORITHM OUTER LOOP: MACS
    algo_outer.optimise = @optimise_macs_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        par_macs.maxnfeval = problem_ebro.dim_d;      % Max function evaluations   100;%
        par_macs.popsize = 1;                              % Population size
        par_macs.rhoini = 1;                                % Initial local hypercube size
        par_macs.F = 1;                                     % F
        par_macs.CR = 0.1;                                  % CR
        par_macs.p_social = 0.5;                            % Ratio between elite and total population
        par_macs.max_arch = 20;                             % Output archive size
        par_macs.coord_ratio = 1;                           % Quota of coordinates examined for each set of individualistic actions
        par_macs.contr_ratio = 0.5;                         % Contraction ratio
        par_macs.draw_flag = 0;                             % Print itarations status  
        par_macs.cp = 0;                                    % constraint flag
        par_macs.MBHflag = 0;                               % number of MBH steps
    algo_outer.par = par_macs;
end



%% ALGORITHM INNER LOOP: MPAIDEA
algo_inner.optimise = @optimise_mpaidea_wrapper;          % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)

        par_mpaidea.nFeValMax = 15000;                          % number of function evaluations for IDEA     100;%
        par_mpaidea.n_populations = 1;                         % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(5,problem_ebro.dim_u);      % number of agents in one population        
        par_mpaidea.max_LR = 20;                               % max.number of local restarts
        par_mpaidea.F = 0.75;                                    % F
        par_mpaidea.CR = 0.8;                                   % CR
        par_mpaidea.parallel = 0;                              % Parallel evaluation of individuals of the DE?        
        par_mpaidea.plots = 0;                                 % Display plots during run?
        par_mpaidea.text = 0;                                  % Display text during run?
        par_mpaidea.fmincon_set = 'sqp';                       % choose algorithm for fmincon: 'sqp', 'interior-point', 
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
        
algo_inner.par = par_mpaidea;






%% ALGORITHM DECOMPOSITION LOOP: MPAIDEA
algo_decomposition.optimise = @optimise_mpaidea_wrapper;          % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)

        par_mpaidea.nFeValMax = 1000*problem_ebro.dim_d;     % number of function evaluations for IDEA     100;%
        par_mpaidea.n_populations = 1;                         % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(5,problem_ebro.dim_u);    % number of agents in one population        
        par_mpaidea.max_LR = 20;                               % max.number of local restarts
        par_mpaidea.F = [];                                    % F
        par_mpaidea.CR = [];                                   % CR
        par_mpaidea.parallel = 0;                              % Parallel evaluation of individuals of the DE?        
        par_mpaidea.plots = 0;                                 % Display plots during run?
        par_mpaidea.text = 0;                                  % Display text during run?
        par_mpaidea.fmincon_set = 'sqp';                       % choose algorithm for fmincon: 'sqp', 'interior-point', 
        %------------------------------------------------------------------
        % default parameters
        %------------------------------------------------------------------
        par_mpaidea.population=[];                              % initial population, same for all execution. Leave empty to randomize.
        par_mpaidea.DE_strategy = 1;                            % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
        par_mpaidea.prob_DE_strategy = 0.5;                     % probability or not using DE/Rand
        par_mpaidea.delta_local = 0.1;                          % dimension of the bubble for the local restart of population.
        par_mpaidea.delta_global = 0.1;                         % characteristic dimension for the global restart of the population
        par_mpaidea.rho = 0.25;                                 % contraction threshold for the population
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
        
algo_decomposition.par = par_mpaidea;
    




%% ------------------------------------------------------------------------
% plot setting
% -------------------------------------------------------------------------
flag = 0;

% algo inner
algo_inner.save_pop_DE = flag;
algo_inner.save_local_search = flag;
algo_inner.save_pop_LR = 0;
algo_inner.save_pop_GR = 0;


algo_inner.save_archive  = flag;

 
% algo outer
algo_outer.save_pop_DE = flag;
algo_outer.save_local_search = flag;
algo_outer.save_pop_LR = 0;
algo_outer.save_pop_GR = 0;
 
% algo constraint
algo_inner.con.save_pop_DE = flag;
algo_inner.con.save_local_search = flag;
algo_inner.con.save_pop_LR = 0;
algo_inner.con.save_pop_GR = 0;

%%-------------------------------------------------------------------------




return