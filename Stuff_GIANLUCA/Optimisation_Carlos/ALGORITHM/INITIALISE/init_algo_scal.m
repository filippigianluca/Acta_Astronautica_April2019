function [algo_minmax, algo_outer, algo_inner] = init_algo_scal(minmax_problem);

%% META ALGORITHM: MACSMINMAX
algo_minmax.optimise = @optimise_scal;            % algorithm in the form [x,f,exitflag,output] = algo(problem,algo_outer,algo_inner,par_minmax)
% parameters
    if isfield(minmax_problem,'maxnfeval')
        par_minmax.maxnfeval = minmax_problem.maxnfeval;% TOTAL function evaluation limit
    else
        error('max nfeval not supplied')
    end
    par_minmax.n_d0 = 4;
    par_minmax.n_subproblems = 8;
    par_minmax.local_search_flags.validation = true;    % Run local search (true) or function evaluation (false) in validation
    par_minmax.local_search_flags.inner = false;        % Idem inner (SO) problem
    par_minmax.local_search_flags.outer = false;       % Idem outer (MO) problem
    % par_minmax.tconv = 1;                               % tconv
    % par_minmax.tspr = 1;                                % tspr
algo_minmax.par_minmax = par_minmax;

% % WITH SO
%     %% ALGORITHM OUTER LOOP: MPAIDEA
%     algo_outer.optimise = @optimise_mpaidea_wrapper;          % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
%         par_mpaidea.nFeValMax = 250*minmax_problem.dim_d;     % number of function evaluations for IDEA
%         par_mpaidea.n_populations = 1;                        % number of populations, if no adaptive behaviour should set to 1
%         par_mpaidea.n_agents = max(5,minmax_problem.dim_d);   % number of agents in one population
%         par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
%         par_mpaidea.max_LR = 10;                              % max.number of local restarts
%         par_mpaidea.DE_strategy = 1;                          % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
%         par_mpaidea.prob_DE_strategy = 0.5;                   % probability or not using DE/Rand
%         par_mpaidea.delta_local = 0.1;                        % dimension of the bubble for the local restart of population.
%         par_mpaidea.delta_global = 0.1;                       % characteristic dimension for the global restart of the population
%         par_mpaidea.rho = 0.20;                               % contraction threshold for the population
%         par_mpaidea.dd_CRF = 3;                               % parameter for the adaptation of CRF.
%         par_mpaidea.F = 0.9;                                  % F
%         par_mpaidea.CR = 0.2;                                 % CR
%         par_mpaidea.text = 0;                                 % verbosity
%         par_mpaidea.record = [1];                             % Fraction(s) of nfevalmax at which the result is saved
%     algo_outer.par = par_mpaidea;


% WITH AUSA
    %% ALGORITHM OUTER LOOP with archive cross-checks
    algo_outer.optimise_au = @optimise_mpaidea_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        par_mpaidea = struct;
        par_mpaidea.nFeValMax = 100*minmax_problem.dim_d;     % number of function evaluations for IDEA
        par_mpaidea.n_populations = 1;                        % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(5,minmax_problem.dim_d);   % number of agents in one population
        par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
        par_mpaidea.max_LR = 10;                              % max.number of local restarts
        par_mpaidea.DE_strategy = 1;                          % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
        par_mpaidea.prob_DE_strategy = 0.5;                   % probability or not using DE/Rand
        par_mpaidea.delta_local = 0.1;                        % dimension of the bubble for the local restart of population.
        par_mpaidea.delta_global = 0.1;                       % characteristic dimension for the global restart of the population
        par_mpaidea.rho = 0.20;                               % contraction threshold for the population
        par_mpaidea.dd_CRF = 3;                               % parameter for the adaptation of CRF.
        par_mpaidea.F = 0.9;                                  % F
        par_mpaidea.CR = 0.2;                                 % CR
        par_mpaidea.text = 0;                                 % verbosity
        par_mpaidea.record = [1];                             % Fraction(s) of nfevalmax at which the result is saved
        par_mpaidea.use = true;
        if(par_mpaidea.nFeValMax == 0)
            par_mpaidea.use = false;
        end
    algo_outer.par_au = par_mpaidea;

    %% ALGORITHM OUTER LOOP with surrogate
    algo_outer.optimise_sa = @optimise_mpaidea_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        %parameters{
        par_mpaidea = struct;
        par_mpaidea.nFeValMax = 1e3*minmax_problem.dim_d;     % number of function evaluations for IDEA
        par_mpaidea.n_populations = 1;                        % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(5,minmax_problem.dim_d);   % number of agents in one population
        par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
        par_mpaidea.max_LR = 10;                              % max.number of local restarts
        par_mpaidea.DE_strategy = 1;                          % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
        par_mpaidea.prob_DE_strategy = 0.5;                   % probability or not using DE/Rand
        par_mpaidea.delta_local = 0.1;                        % dimension of the bubble for the local restart of population.
        par_mpaidea.delta_global = 0.1;                       % characteristic dimension for the global restart of the population
        par_mpaidea.rho = 0.20;                               % contraction threshold for the population
        par_mpaidea.dd_CRF = 3;                               % parameter for the adaptation of CRF.
        par_mpaidea.F = 0.9;                                  % F
        par_mpaidea.CR = 0.2;                                 % CR
        par_mpaidea.text = 0;                                 % verbosity
        par_mpaidea.record = [1];                             % Fraction(s) of nfevalmax at which the result is saved
        par_mpaidea.use = true;
        if(par_mpaidea.nFeValMax == 0)
            par_mpaidea.use = false;
        end
    algo_outer.par_sa = par_mpaidea;

%% ALGORITHM INNER LOOP: MPAIDEA
algo_inner.optimise = @optimise_mpaidea_wrapper;          % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    par_mpaidea.nFeValMax = 600*minmax_problem.dim_u;   % number of function evaluations for IDEA
    par_mpaidea.n_populations = 1;                        % number of populations, if no adaptive behaviour should set to 1
    par_mpaidea.n_agents = max(5,minmax_problem.dim_u);   % number of agents in one population
    par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
    par_mpaidea.max_LR = 10;                               % max.number of local restarts
    par_mpaidea.DE_strategy = 1;                          % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
    par_mpaidea.prob_DE_strategy = 0.5;                   % probability or not using DE/Rand
    par_mpaidea.delta_local = 0.1;                        % dimension of the bubble for the local restart of population.
    par_mpaidea.delta_global = 0.1;                       % characteristic dimension for the global restart of the population
    par_mpaidea.rho = 0.20;                               % contraction threshold for the population
    par_mpaidea.dd_CRF = 3;                               % parameter for the adaptation of CRF.
    par_mpaidea.F = 0.9;                                  % F
    par_mpaidea.CR = 0.2;                                 % CR
    par_mpaidea.text = 0;                                 % verbosity
    par_mpaidea.record = [1];                        % Fraction(s) of nfevalmax at which the result is saved
algo_inner.par = par_mpaidea;


return