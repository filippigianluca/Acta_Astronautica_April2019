function [algo_inner] = init_algo_relax_bisurr_advanced_checks(minmax_problem);


algo_inner = cell(1,minmax_problem.n_obj);
for obj = 1:minmax_problem.n_obj

    % %% ALGORITHM INNER LOOP: FMINCON's SQP
    % algo_inner{obj}.optimise = @optimise_fmincon_wrapper;                % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    % par_local = struct;
    %     % default parameters --> random choice of u0
    % algo_inner{obj}.par = par_local;

    %% ALGORITHM INNER LOOP: MPAIDEA
    algo_inner{obj}.optimise = @optimise_mpaidea_wrapper;                % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    par_mpaidea = struct;
        par_mpaidea.nFeValMax = 5e4*minmax_problem.dim_u;     % number of function evaluations for IDEA
        par_mpaidea.n_populations = 4;                        % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = 10;%max(4,minmax_problem.dim_u);   % number of agents in one population
        par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
        par_mpaidea.max_LR = 5;                               % max.number of local restarts
        par_mpaidea.DE_strategy = 1;                          % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
        par_mpaidea.prob_DE_strategy = 0.5;                   % probability or not using DE/Rand
        par_mpaidea.delta_local = 0.1;                        % dimension of the bubble for the local restart of population.
        par_mpaidea.delta_global = 0.1;                       % characteristic dimension for the global restart of the population
        par_mpaidea.rho = 0.25;                               % contraction threshold for the population
        par_mpaidea.dd_CRF = 3;                               % parameter for the adaptation of CRF.
        par_mpaidea.F = [];                                  % F
        par_mpaidea.CR = [];                                 % CR
        par_mpaidea.text = 1;                                 % verbosity
        par_mpaidea.record = 0.05:0.05:1;% Fraction(s) of nfevalmax at which the result is saved
    algo_inner{obj}.par = par_mpaidea;


end



return