function [algo_minmax, algo_outer, algo_inner] = init_algo_so_nested_local(minmax_problem);

%% META ALGORITHM: MACSMINMAX
algo_minmax.optimise = @optimise_nested;            % algorithm in the form [x,f,exitflag,output] = algo(problem,algo_outer,algo_inner,par_minmax)
% parameters
    if isfield(minmax_problem,'maxnfeval')
        maxnfeval = minmax_problem.maxnfeval*100;% TOTAL function evaluation limit
        par_minmax.maxnfeval = maxnfeval;
    else
        error('max nfeval not supplied')
    end

    nfeval_outer = [];
    nfeval_inner = [];

    if ~isempty(nfeval_outer) && ~isempty(nfeval_inner)
        warning('Supplied nfevalmax, nfeval_outer and nfeval_inner => will ignore nfevalmax')
    elseif ~isempty(nfeval_inner)
        nfeval_outer = floor(maxnfeval/nfeval_inner);
    elseif ~isempty(nfeval_outer)
        nfeval_inner = floor(maxnfeval/nfeval_outer);
    else
        nfeval_outer = floor(sqrt(maxnfeval/minmax_problem.n_obj));
        nfeval_inner = nfeval_outer;
    end

algo_minmax.par_minmax = par_minmax;

if (minmax_problem.n_obj == 1)
    %% ALGORITHM OUTER LOOP: MPAIDEA
    algo_outer.optimise = @optimise_mpaidea_wrapper;          % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
        par_mpaidea.nFeValMax = minmax_problem.nfouter_nested*minmax_problem.dim_d;% number of function evaluations for IDEA
        par_mpaidea.n_populations = 1;                        % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(4,minmax_problem.dim_d);   % number of agents in one population
        par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
        par_mpaidea.max_LR = 3;                               % max.number of local restarts
        par_mpaidea.DE_strategy = 1;                          % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
        par_mpaidea.prob_DE_strategy = 0.5;                   % probability or not using DE/Rand
        par_mpaidea.delta_local = 0.1;                        % dimension of the bubble for the local restart of population.
        par_mpaidea.delta_global = 0.1;                       % characteristic dimension for the global restart of the population
        par_mpaidea.rho = 0.25;                               % contraction threshold for the population
        par_mpaidea.dd_CRF = 3;                               % parameter for the adaptation of CRF.
        par_mpaidea.F = 0.2;                                  % F
        par_mpaidea.CR = 0.8;                                 % CR
        par_mpaidea.text = 0;                                 % verbosity
        par_mpaidea.record = [1];                             % Fraction(s) of nfevalmax at which the result is saved
    algo_outer.par = par_mpaidea;


else
    %% ALGORITHM OUTER LOOP: MACS
    algo_outer.optimise = @optimise_macs_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        par_macs.maxnfeval = nfeval_outer;                  % Max function evaluations
        par_macs.popsize = 10 + minmax_problem.dim_d;       % Population size
        par_macs.rhoini = 1;                                % Initial local hypercube size
        par_macs.F = 0.9;                                   % F
        par_macs.CR = 0.9;                                  % CR
        par_macs.p_social = 0.2;                            % Ratio between elite and total population
        par_macs.max_arch = 100;                            % Inner archive size
        par_macs.max_arch_out = 40;                          % Output archive size
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
        par_macs.optimal_control = 1;
        par_macs.vars_to_opt = logical([ones(minmax_problem.dim_d,1);zeros(minmax_problem.n_obj*minmax_problem.dim_u,1)]);


    algo_outer.par = par_macs;
end

algo_inner = cell(1,minmax_problem.n_obj);
for obj = 1:minmax_problem.n_obj
    %% ALGORITHM INNER LOOP: FMINCON's SQP
    algo_inner{obj}.optimise = @optimise_fmincon_wrapper;                % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    par_local = struct;
        % default parameters --> random choice of u0
    algo_inner{obj}.par = par_local;

    % %% ALGORITHM INNER LOOP: MPAIDEA
    % algo_inner{obj}.optimise = @optimise_mpaidea_wrapper;                % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    % par_mpaidea = struct;
    %     par_mpaidea.nFeValMax = nfeval_inner;           % number of function evaluations for IDEA
    %     par_mpaidea.n_populations = 1;                              % number of populations, if no adaptive behaviour should be set to 1
    %     par_mpaidea.n_agents = max(4,minmax_problem.dim_u);         % number of agents in one population
    %     par_mpaidea.initial_population=[];                          % initial population, same for all execution. Leave empty to randomize.
    %     par_mpaidea.max_LR = 1;                                     % max.number of local restarts
    %     par_mpaidea.DE_strategy = 1;                                % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
    %     par_mpaidea.prob_DE_strategy = 0.50;                        % probability of not using DE/Rand
    %     par_mpaidea.delta_local = 0.10;                             % dimension of the bubble for the local restart of population.
    %     par_mpaidea.delta_global = 0.10;                            % characteristic dimension for the global restart of the population
    %     par_mpaidea.rho = 0.20;                                     % contraction threshold for the population
    %     par_mpaidea.dd_CRF = 3;                                     % parameter for the adaptation of CRF.
    %     par_mpaidea.F = 0.90;                                       % F
    %     par_mpaidea.CR = 0.20;                                      % CR
    %     par_mpaidea.text = 0;                                       % verbosity
    %     par_mpaidea.record = [1];                                   % Fraction(s) of nfevalmax at which the result is saved
    % algo_inner{obj}.par = par_mpaidea;
    
    % %% ALGORITHM INNER LOOP: IDEA2
    % algo_inner{obj}.optimise = @optimise_idea2_wrapper;            % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    % %parameters{
    % options=[];
    %     options_idea2(1) = nfeval_inner;                     % number of function evaluations for IDEA
    %     % options_idea2(1) = 500;
    %     options_idea2(4) = max(4,minmax_problem.dim_u);       % number of agents
    %     options_idea2(31) = 20;                               % number of local restarts
    %     options_idea2(32) = 2;                                % communication heuristics (see com_red.m)
    %     options_idea2(22) = 0.1;                              % delta_c: crowding factor for global restart
    %     options_idea2(27) = 0.1;                              % tolconv: local convergence
    %     options_idea2(15) = -1;                               % verbousity
    %     options_idea2(28) = 1;                                % F
    %     options_idea2(29) = 0.1;                              % CR
    %     % options_idea2(30) = minmax_problem.n_obj;             % number of objectives (not used)
    %     options_idea2(33) = 0.2;                              % initial probability of running the subproblem (not used)
        
    % algo_inner{obj}.par.options = options_idea2;
    % algo_inner{obj}.par.initial_population = [];
end

return