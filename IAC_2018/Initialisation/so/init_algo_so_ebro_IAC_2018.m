% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [algo_minmax, algo_outer, algo_inner, algo_constraint] = init_algo_so_ebro_IAC_2018(minmax_problem)
%
% [algo_minmax, algo_outer, algo_inner] = init_algo_so(minmax_problem)
%
% define settings for the optimisers in the minmax algorithm: 
% optimisaer for: inner loop, outer loop, constraint maximisation.
% choose to set parameters for MPAIDEA or for fmincon




%======================================================================================================================================================
%% META ALGORITHM: MACSMINMAX
algo_minmax.optimise = @optimise_minmax_so;            % algorithm in the form [x,f,exitflag,output] = algo(problem,algo_outer,algo_inner,par_minmax)
% parameters
    if isfield(minmax_problem,'maxnfeval')
        par_minmax.maxnfeval = minmax_problem.maxnfeval;% TOTAL function evaluation limit
    else
        par_minmax.maxnfeval = 2e5;                     %error('max nfeval not supplied') % 100;%
    end
    par_minmax.n_d0 = 1;
    
    par_minmax.local_search_flags.validation = false;    % Run local search (true) or function evaluation (false) in validation
    par_minmax.local_search_flags.inner = false;        % Idem inner (SO) problem
    par_minmax.local_search_flags.outer = false;       % Idem outer (MO) problem

algo_minmax.par_minmax = par_minmax;
%======================================================================================================================================================





%======================================================================================================================================================
    %% ALGORITHM OUTER LOOP: MPAIDEA
    algo_outer.optimise = @optimise_mpaidea_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
     
        par_mpaidea.nFeValMax = 5000;                          % number of function evaluations for IDEA     100;%
        par_mpaidea.n_populations = 1;                         % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(5,minmax_problem.dim_d);    % number of agents in one population        
        par_mpaidea.max_LR = 10;                               % max.number of local restarts
        par_mpaidea.F  = [];                                   % F
        par_mpaidea.CR = [];                                   % CR
        par_mpaidea.parallel = 0;                              % Parallel evaluation of individuals of the DE?        
        par_mpaidea.plots = 0;                                 % Display plots during run?
        par_mpaidea.text = 0;                                  % Display text during run?
        par_mpaidea.fmincon_set = 'sqp';                       % choose algorithm for fmincon: 'sqp', 'interior-point', 
        par_mpaidea.fitnessfcn.obj_constr = minmax_problem.obj_constr;                  % Flag to 1 if objective and constraints are in the same function
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
        par_mpaidea.record = 0:0.01:1;                               % Fraction(s) of nfevalmax at which the result is saved       
        par_mpaidea.save_pop_DE = 0;                            % Save results of DE to file?
        par_mpaidea.name_save_pop_DE = 'pop_DE_';               % If yes, choose prefix of name for files:
        par_mpaidea.save_local_search = 0;                      % Save results of local search to file?
        par_mpaidea.name_save_local_search = 'minima_fmincon_'; % If yes, choose prefix for name for file:
        par_mpaidea.save_pop_LR = 0;                            % Save populations at local restart (each one saved on a different file)?
        par_mpaidea.name_save_pop_LR = 'pop_LR_';               % If yes, choose prefix of name for files:
        par_mpaidea.save_pop_GR = 0;                            % Save populations at global restart (each one saved on a different file)?
        par_mpaidea.name_save_pop_GR = 'pop_GR_';               % If yes, choose prefix of name for files
        
        par_mpaidea.fitnessfcn.weighted   = 0;                  % How to handle constraints: set to 1 for weighted constraints with fixed weights, or to 0 for penalty with no weights
        par_mpaidea.fitnessfcn.ceq_eps    = 1e-6;               % If the constraints are handled without weights, then define a tolerance for the violation of the equality constraints
        par_mpaidea.fitnessfcn.w_ceq      = 1000;               % Weights for penalty
        par_mpaidea.fitnessfcn.w_c        = 100;
        
  algo_outer.par = par_mpaidea;


% %% ALGORITHM OUTER LOOP: IDEA2
% algo_outer.optimise = @optimise_idea2_wrapper;            % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
% %parameters{
%     options_idea2(1) = 200*minmax_problem.dim_u;          % number of function evaluations for IDEA
%     options_idea2(4) = 5;                                 % number of agents
%     options_idea2(31) = 100;                              % number of local restarts
%     options_idea2(32) = 2;                                % communication heuristics (see com_red.m)
%     options_idea2(22) = 0.1;                              % delta_c: crowding factor for global restart
%     options_idea2(27) = 0.25;                             % tolconv: local convergence
%     options_idea2(15) = -1;                               % verbousity
%     options_idea2(28) = 0.2;                                % F
%     options_idea2(29) = 0.8;                                % CR
%     % options_idea2(30) = minmax_problem.n_obj;             % number of objectives (not used)
%     options_idea2(33) = 0.2;                              % initial probability of running the subproblem (not used)

% algo_outer.par.options = options_idea2;
% algo_outer.par.initial_population = []; 
    
       
    
% %% ALGORITHM OUTER LOOP: fmincon
%     algo_outer.optimise = @optimise_fmincon_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)  
%     
%         par_fmincon.opt        = 'sqp';    
%         par_fmincon.initial_population = [];
%         par_fmincon.nFeValMax = 1000*minmax_problem.dim_d;
%         
%         par_fmincon.fitnessfcn.obj_constr = minmax_problem.obj_constr;  
%   algo_outer.par = par_fmincon;
%======================================================================================================================================================
  
  
  









%======================================================================================================================================================
%% ALGORITHM INNER LOOP: MPAIDEA
algo_inner.optimise = @optimise_mpaidea_wrapper;          % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)

        par_mpaidea.nFeValMax = 5000;     % number of function evaluations for IDEA     100;%
        par_mpaidea.n_populations = 1;                         % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(5,minmax_problem.dim_u);    % number of agents in one population        
        par_mpaidea.max_LR = 10;                               % max.number of local restarts
        par_mpaidea.F  = [];                                   % F
        par_mpaidea.CR = [];                                   % CR
        par_mpaidea.parallel = 0;                              % Parallel evaluation of individuals of the DE?        
        par_mpaidea.plots = 0;                                 % Display plots during run?
        par_mpaidea.text = 0;                                  % Display text during run?
        par_mpaidea.fmincon_set = 'sqp';                       % choose algorithm for fmincon: 'sqp', 'interior-point', 
        par_mpaidea.fitnessfcn.obj_constr = minmax_problem.obj_constr;                  % Flag to 1 if objective and constraints are in the same function
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
        par_mpaidea.record = 0:0.01:1;                               % Fraction(s) of nfevalmax at which the result is saved       
        par_mpaidea.save_pop_DE = 0;                            % Save results of DE to file?
        par_mpaidea.name_save_pop_DE = 'pop_DE_';               % If yes, choose prefix of name for files:
        par_mpaidea.save_local_search = 0;                      % Save results of local search to file?
        par_mpaidea.name_save_local_search = 'minima_fmincon_'; % If yes, choose prefix for name for file:
        par_mpaidea.save_pop_LR = 0;                            % Save populations at local restart (each one saved on a different file)?
        par_mpaidea.name_save_pop_LR = 'pop_LR_';               % If yes, choose prefix of name for files:
        par_mpaidea.save_pop_GR = 0;                            % Save populations at global restart (each one saved on a different file)?
        par_mpaidea.name_save_pop_GR = 'pop_GR_';               % If yes, choose prefix of name for files
        
        par_mpaidea.fitnessfcn.weighted   = 0;                  % How to handle constraints: set to 1 for weighted constraints with fixed weights, or to 0 for penalty with no weights
        par_mpaidea.fitnessfcn.ceq_eps    = 1e-6;               % If the constraints are handled without weights, then define a tolerance for the violation of the equality constraints
        par_mpaidea.fitnessfcn.w_ceq      = 1000;               % Weights for penalty
        par_mpaidea.fitnessfcn.w_c        = 100;
        
algo_inner.par = par_mpaidea;




% %% ALGORITHM INNER LOOP: IDEA2
% algo_inner.optimise = @optimise_idea2_wrapper;            % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
% %parameters{
%     options_idea2(1) = 200*minmax_problem.dim_u;          % number of function evaluations for IDEA
%     options_idea2(4) = 5;                                 % number of agents
%     options_idea2(31) = 100;                              % number of local restarts
%     options_idea2(32) = 2;                                % communication heuristics (see com_red.m)
%     options_idea2(22) = 0.1;                              % delta_c: crowding factor for global restart
%     options_idea2(27) = 0.25;                             % tolconv: local convergence
%     options_idea2(15) = -1;                               % verbousity
%     options_idea2(28) = 0.2;                                % F
%     options_idea2(29) = 0.8;                                % CR
%     % options_idea2(30) = minmax_problem.n_obj;             % number of objectives (not used)
%     options_idea2(33) = 0.2;                              % initial probability of running the subproblem (not used)

% algo_inner.par.options = options_idea2;
% algo_inner.par.initial_population = [];




% %% ALGORITHM INNER LOOP: fmincon
%     algo_inner.optimise = @optimise_fmincon_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)  
%     
%         par_fmincon.opt        = 'sqp';    
%         par_fmincon.initial_population = [];
%         par_fmincon.nFeValMax = 1000*minmax_problem.dim_u;
%         
%         par_fmincon.fitnessfcn.obj_constr = minmax_problem.obj_constr; 
%   algo_inner.par = par_fmincon;
%======================================================================================================================================================


  
  
  
  
  
  
  
  
  
  
%======================================================================================================================================================  
%% ALGORITHM INNER LOOP (CONSTRAINT): MPAIDEA
algo_constraint.optimise = @optimise_mpaidea_wrapper;          % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)

        par_mpaidea.nFeValMax = 5000;     % number of function evaluations for IDEA     100;%
        par_mpaidea.n_populations = 1;                         % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(5,minmax_problem.dim_u);    % number of agents in one population        
        par_mpaidea.max_LR = 10;                               % max.number of local restarts
        par_mpaidea.F  = [];                                   % F
        par_mpaidea.CR = [];                                   % CR
        par_mpaidea.parallel = 0;                              % Parallel evaluation of individuals of the DE?        
        par_mpaidea.plots = 0;                                 % Display plots during run?
        par_mpaidea.text = 0;                                  % Display text during run?
        par_mpaidea.fmincon_set = 'sqp';                       % choose algorithm for fmincon: 'sqp', 'interior-point', 
        par_mpaidea.fitnessfcn.obj_constr = minmax_problem.obj_constr;                  % Flag to 1 if objective and constraints are in the same function
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
        par_mpaidea.record = 0:0.01:1;                               % Fraction(s) of nfevalmax at which the result is saved       
        par_mpaidea.save_pop_DE = 0;                            % Save results of DE to file?
        par_mpaidea.name_save_pop_DE = 'pop_DE_';               % If yes, choose prefix of name for files:
        par_mpaidea.save_local_search = 0;                      % Save results of local search to file?
        par_mpaidea.name_save_local_search = 'minima_fmincon_'; % If yes, choose prefix for name for file:
        par_mpaidea.save_pop_LR = 0;                            % Save populations at local restart (each one saved on a different file)?
        par_mpaidea.name_save_pop_LR = 'pop_LR_';               % If yes, choose prefix of name for files:
        par_mpaidea.save_pop_GR = 0;                            % Save populations at global restart (each one saved on a different file)?
        par_mpaidea.name_save_pop_GR = 'pop_GR_';               % If yes, choose prefix of name for files
        
        par_mpaidea.fitnessfcn.weighted   = 0;                  % How to handle constraints: set to 1 for weighted constraints with fixed weights, or to 0 for penalty with no weights
        par_mpaidea.fitnessfcn.ceq_eps    = 1e-6;               % If the constraints are handled without weights, then define a tolerance for the violation of the equality constraints
        par_mpaidea.fitnessfcn.w_ceq      = 1000;               % Weights for penalty
        par_mpaidea.fitnessfcn.w_c        = 100;
        par_mpaidea.fitnessfcn.obj_constr = 0; 
        par_mpaidea.fitnessfcn.constr     = [];
        
        
        
algo_constraint.par = par_mpaidea;






% %% ALGORITHM INNER LOOP (CONSTRAINT): fmincon
%     algo_constraint.optimise = @optimise_fmincon_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)  
%     
%         par_fmincon.opt        = 'sqp';    
%         par_fmincon.initial_population = [];
%         par_fmincon.nFeValMax = 1000*minmax_problem.dim_u;
%         
%         par_fmincon.fitnessfcn.obj_constr = minmax_problem.obj_constr; 
%   algo_constraint.par = par_fmincon;
%======================================================================================================================================================



  
  
  
  
  
  
  
  
  
  
  

%%
%%-------------------------------------------------------------------------

% PLOT

% algo inner
 
algo_inner.save_pop_DE = 0;
algo_inner.save_local_search = 0;
algo_inner.save_pop_LR = 1;
algo_inner.save_pop_GR = 1;
 
% algo outer
 
algo_outer.save_pop_DE = 0;
algo_outer.save_local_search = 0;
algo_outer.save_pop_LR = 0;
algo_outer.save_pop_GR = 0;
 
% algo constraint
 
algo_inner.con.save_pop_DE = 0;
algo_inner.con.save_local_search = 0;
algo_inner.con.save_pop_LR = 0;
algo_inner.con.save_pop_GR = 0;



% -------------------------------------------------------------------------
% plot
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