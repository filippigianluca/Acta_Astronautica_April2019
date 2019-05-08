function [algo_inner, algo_outer, algo_minmax, algo_decomposition] = def_optimiser_IAC_2018_so(problem)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MINMAX:
%
% MP_AIDEA Parameters


% META LOOP

%--------------------------------------------------------------------------
% algorithm in the form [x,f,exitflag,output] = algo(problem,algo_outer,algo_inner,par_minmax)
%--------------------------------------------------------------------------
algo_minmax.optimise = @optimise_so;

%--------------------------------------------------------------------------
% maximum number of function evaluation
%--------------------------------------------------------------------------
par_minmax.maxnfeval = 1e6;

%--------------------------------------------------------------------------
% number of initial design vectors
%--------------------------------------------------------------------------
par_minmax.n_d0 = 1;

%--------------------------------------------------------------------------
% Run local search (true) or function evaluation (false) in validation
%--------------------------------------------------------------------------
par_minmax.local_search_flags.validation = false;
par_minmax.local_search_flags.inner = false;       % Idem inner (SO) problem
par_minmax.local_search_flags.outer = false;       % Idem outer (MO) problem

algo_minmax.par_minmax = par_minmax;





% INNER LOOP (maximisation over u)

%--------------------------------------------------------------------------
% algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
%--------------------------------------------------------------------------
algo_inner.optimise = @optimise_mpaidea_wrapper;

%--------------------------------------------------------------------------
% maximum number of function evaluation for the inner loop
%--------------------------------------------------------------------------
par_mpaidea.nFeValMax = 2000;

%--------------------------------------------------------------------------
% number of populations, if no adaptive behaviour should set to 1
%--------------------------------------------------------------------------
par_mpaidea.n_populations = 1;

%--------------------------------------------------------------------------
% number of agents in one population
%--------------------------------------------------------------------------
par_mpaidea.n_agents = max(5,problem.dim_u);

% -------------------------------------------------------------------------
% Maximum number of local restart before global restart (for cases when
% only one population is considered and no adaptation of delta_local and
% local/global restart is performed)
% par_mpaidea.max_LR = []; --> adaptation
% -------------------------------------------------------------------------
par_mpaidea.max_LR = 20;

% -------------------------------------------------------------------------
% other parameters: default values
% -------------------------------------------------------------------------
par_mpaidea.parallel = 0;                               % Parallel evaluation of individuals of the DE?
par_mpaidea.plots = 0;                                  % Display plots during run?
par_mpaidea.text = 0;                                   % Display text during run?
par_mpaidea.save_pop_DE = 0;                            % Save results of DE to file?
par_mpaidea.name_save_pop_DE = 'pop_DE_spacecraft_';               % If yes, choose prefix of name for files:
par_mpaidea.save_local_search = 0;                      % Save results of local search to file?
par_mpaidea.name_save_local_search = 'minima_fmincon_spacecraft_'; % If yes, choose prefix for name for file:
par_mpaidea.save_pop_LR = 0;                            % Save populations at local restart (each one saved on a different file)?
name_save_pop_LR = 'pop_LR_spacecraft_';                           % If yes, choose prefix of name for files:
par_mpaidea.save_pop_GR = 0;                            % Save populations at global restart (each one saved on a different file)?
name_save_pop_GR = 'pop_GR_spacecraft_';                           % If yes, choose prefix of name for files:

par_mpaidea.population=[];                     % initial population, same for all execution. Leave empty to randomize.
par_mpaidea.DE_strategy = 1;                   % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
par_mpaidea.prob_DE_strategy = 0.5;            % probability or not using DE/Rand
par_mpaidea.delta_local = 0.1;                 % dimension of the bubble for the local restart of population.
par_mpaidea.delta_global = 0.1;                % characteristic dimension for the global restart of the population
par_mpaidea.rho = 0.25;                        % contraction threshold for the population
par_mpaidea.dd_CRF = 3;                        % parameter for the adaptation of CRF.
par_mpaidea.F = [];                           % F
par_mpaidea.CR = [];                          % CR
par_mpaidea.text = 0;                          % verbosity
par_mpaidea.record = [1];                      % Fraction(s) of nfevalmax at which the result is saved


% -------------------------------------------------------------------------
% set for fmincon: local search before local/globa restart
if isempty(problem.constraints{1})
    par_mpaidea.fmincon_set = 'sqp';
else
    par_mpaidea.fmincon_set = 'interior-point';
end
% -------------------------------------------------------------------------

algo_inner.par = par_mpaidea;



% OUTER LOOP (minimisation over d)

%--------------------------------------------------------------------------
% algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
%--------------------------------------------------------------------------
algo_outer.optimise = @optimise_mpaidea_wrapper;

%--------------------------------------------------------------------------
% maximum number of function evaluation for the outer loop
%--------------------------------------------------------------------------
par_mpaidea.nFeValMax = 2000;

%--------------------------------------------------------------------------
% number of populations, if no adaptive behaviour should set to 1
%--------------------------------------------------------------------------
par_mpaidea.n_populations = 1;

%--------------------------------------------------------------------------
% number of agents in one population
%--------------------------------------------------------------------------
par_mpaidea.n_agents = max(5,problem.dim_d);

% -------------------------------------------------------------------------
% Maximum number of local restart before global restart (for cases when
% only one population is considered and no adaptation of delta_local and
% local/global restart is performed)
% par_mpaidea.max_LR = []; --> adaptation
% -------------------------------------------------------------------------
par_mpaidea.max_LR = 20;

% -------------------------------------------------------------------------
% other parameters: default values
% -------------------------------------------------------------------------
par_mpaidea.parallel = 0;                               % Parallel evaluation of individuals of the DE?
par_mpaidea.plots = 0;                                  % Display plots during run?
par_mpaidea.text = 0;                                   % Display text during run?
par_mpaidea.save_pop_DE = 0;                            % Save results of DE to file?
par_mpaidea.name_save_pop_DE = 'pop_DE_spacecraft_';               % If yes, choose prefix of name for files:
par_mpaidea.save_local_search = 0;                      % Save results of local search to file?
par_mpaidea.name_save_local_search = 'minima_fmincon_spacecraft_'; % If yes, choose prefix for name for file:
par_mpaidea.save_pop_LR = 0;                            % Save populations at local restart (each one saved on a different file)?
name_save_pop_LR = 'pop_LR_spacecraft_';                           % If yes, choose prefix of name for files:
par_mpaidea.save_pop_GR = 0;                            % Save populations at global restart (each one saved on a different file)?
name_save_pop_GR = 'pop_GR_spacecraft_';                           % If yes, choose prefix of name for files:

par_mpaidea.population=[];                % initial population, same for all execution. Leave empty to randomize.
par_mpaidea.DE_strategy = 1;              % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
par_mpaidea.prob_DE_strategy = 0.5;       % probability or not using DE/Rand
par_mpaidea.delta_local = 0.1;            % dimension of the bubble for the local restart of population.
par_mpaidea.delta_global = 0.1;           % characteristic dimension for the global restart of the population
par_mpaidea.rho = 0.25;                   % contraction threshold for the population
par_mpaidea.dd_CRF = 3;                   % parameter for the adaptation of CRF.
par_mpaidea.F = [];                      % F
par_mpaidea.CR = [];                     % CR
par_mpaidea.text = 0;                     % verbosity
par_mpaidea.record = [1];                 % Fraction(s) of nfevalmax at which the result is saved

% -------------------------------------------------------------------------
% set for fmincon: local search before local/globa restart
if isempty(problem.constraints{1})
    par_mpaidea.fmincon_set = 'sqp';
else
    par_mpaidea.fmincon_set = 'interior-point';
end
% -------------------------------------------------------------------------

algo_outer.par = par_mpaidea;









% DECOMPOSITION ALGORITHM

%--------------------------------------------------------------------------
% algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
%--------------------------------------------------------------------------
algo_decomposition.optimise = @optimise_mpaidea_wrapper;

%--------------------------------------------------------------------------
% maximum number of function evaluation for the inner loop
%--------------------------------------------------------------------------
par_mpaidea.nFeValMax = 5;

%--------------------------------------------------------------------------
% number of populations, if no adaptive behaviour should set to 1
%--------------------------------------------------------------------------
par_mpaidea.n_populations = 1;

%--------------------------------------------------------------------------
% number of agents in one population
%--------------------------------------------------------------------------
par_mpaidea.n_agents = max(5,problem.dim_u);

% -------------------------------------------------------------------------
% Maximum number of local restart before global restart (for cases when
% only one population is considered and no adaptation of delta_local and
% local/global restart is performed)
% par_mpaidea.max_LR = []; --> adaptation
% -------------------------------------------------------------------------
par_mpaidea.max_LR = 20;

% -------------------------------------------------------------------------
% other parameters: default values
% -------------------------------------------------------------------------
par_mpaidea.parallel = 0;                                   % Parallel evaluation of individuals of the DE?
par_mpaidea.plots = 0;                                  % Display plots during run?
par_mpaidea.text = 0;                                   % Display text during run?
par_mpaidea.save_pop_DE = 0;                            % Save results of DE to file?
par_mpaidea.name_save_pop_DE = 'pop_DE_';               % If yes, choose prefix of name for files:
par_mpaidea.save_local_search = 0;                      % Save results of local search to file?
par_mpaidea.name_save_local_search = 'minima_fmincon_'; % If yes, choose prefix for name for file:
par_mpaidea.save_pop_LR = 0;                            % Save populations at local restart (each one saved on a different file)?
name_save_pop_LR = 'pop_LR_';                           % If yes, choose prefix of name for files:
par_mpaidea.save_pop_GR = 0;                            % Save populations at global restart (each one saved on a different file)?
name_save_pop_GR = 'pop_GR_';                           % If yes, choose prefix of name for files:

par_mpaidea.population=[];                     % initial population, same for all execution. Leave empty to randomize.
par_mpaidea.DE_strategy = 1;                   % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
par_mpaidea.prob_DE_strategy = 0.5;            % probability or not using DE/Rand
par_mpaidea.delta_local = 0.1;                 % dimension of the bubble for the local restart of population.
par_mpaidea.delta_global = 0.1;                % characteristic dimension for the global restart of the population
par_mpaidea.rho = 0.25;                        % contraction threshold for the population
par_mpaidea.dd_CRF = 3;                        % parameter for the adaptation of CRF.
par_mpaidea.F = [];                           % F
par_mpaidea.CR = [];                          % CR
par_mpaidea.text = 0;                          % verbosity
par_mpaidea.record = [1];                      % Fraction(s) of nfevalmax at which the result is saved
% -------------------------------------------------------------------------
% set for fmincon: local search before local/globa restart
if isempty(problem.constraints{1})
    par_mpaidea.fmincon_set = 'sqp';
else
    par_mpaidea.fmincon_set = 'interior-point';
end
% -------------------------------------------------------------------------

algo_decomposition.par = par_mpaidea;




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







end