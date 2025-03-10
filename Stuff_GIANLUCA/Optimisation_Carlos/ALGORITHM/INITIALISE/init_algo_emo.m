function [algo_minmax, algo_outer, algo_inner] = init_algo_emo(minmax_problem, varargin);

%% META ALGORITHM: MACSMINMAX
algo_minmax.optimise = @optimise_minmarek;            % algorithm in the form [x,f,exitflag,output] = algo(problem,algo_outer,algo_inner,par_minmax)
% parameters
    if isfield(minmax_problem,'maxnfeval')
        par_minmax.maxnfeval = minmax_problem.maxnfeval/10;% function evaluation limit (does not account for ls in validation)
    else
        error('max nfeval not supplied')
    end

    par_minmax.local_search_flags.validation = false;
    par_minmax.keep_d_record = true;
    par_minmax.tol_conv = 1e-3; %Tolerance to convergence, only used in S.O. case
    par_minmax.indicator_d_max = 1e-3; % Threshold on the indicator (EI, PI, etc.) for outer loop
    par_minmax.indicator_u_max = 1e-3; % Threshold on the indicator (EI, PI, etc.) for inner loop
    par_minmax.iter_d_max = 20; %20*minmax_problem.dim_d; % Max. iterations for outer loop in an iteration of MinMaReK
    par_minmax.iter_u_max = 20; %20*minmax_problem.dim_u; % Max. iterations for inner loop in an iteration of MinMaReK

    
    % initialisation of the surrogate parameters
        if minmax_problem.n_obj == 1
            par_minmax.indicator_d = str2func('kriging_EI2');
        elseif minmax_problem.n_obj == 2
            par_minmax.indicator_d = str2func('kriging_EIaug');
        else
            error('We do not have kriging indicators for many-objective (yet)' );
        end
        
        surrogate = cell(1,minmax_problem.n_obj);
        for obj = 1:minmax_problem.n_obj
            map_u_info = get_map_info(minmax_problem.lb_u{obj}, minmax_problem.ub_u{obj});

            num_fe = 1;
            for dim = 1:minmax_problem.dim_u
                num_fe = num_fe * map_u_info.n_int{dim};
            end

            surrogate{obj} = struct;
            surrogate{obj}.dim_d = minmax_problem.dim_d;
            surrogate{obj}.dim_u = minmax_problem.dim_u;
            surrogate{obj}.method = 'kriging';
            surrogate{obj}.corrfun = @corrgauss;
            surrogate{obj}.regrfun = @regpoly0;
            surrogate{obj}.training = str2func([lower(surrogate{obj}.method) '_training']);
            surrogate{obj}.predictor = str2func([lower(surrogate{obj}.method) '_predictor']);
            surrogate{obj}.indicator_u = str2func([lower(surrogate{obj}.method) '_EI']);

            surrogate{obj}.add = @surrogate_add;
            surrogate{obj}.update = @surrogate_update;
            if (par_minmax.multiple_surrogates)
                surrogate{obj}.model = cell(1,num_fe);
                surrogate{obj}.changed_FE = zeros(1,num_fe);
                surrogate{obj}.x_doe = cell(1,num_fe);
                surrogate{obj}.f_doe = cell(1,num_fe);
                surrogate{obj}.map_info = map_u_info;
                
                surrogate{obj}.find = @surrogate_find_multiple;
                surrogate{obj}.num_FE = num_fe;
            else
                surrogate{obj}.model = cell(1,1);
                surrogate{obj}.changed_FE = ones(1,1);
                surrogate{obj}.x_doe = cell(1,1);
                surrogate{obj}.f_doe = cell(1,1);
                surrogate{obj}.map_info = map_u_info;
                
                surrogate{obj}.find = @surrogate_find_one;
                surrogate{obj}.num_FE = 1;    
            end
        end
    par_minmax.surrogate = surrogate;

algo_minmax.par_minmax = par_minmax;


%% ALGORITHM OUTER LOOP: IDEA2
algo_outer.optimise = @optimise_idea2_wrapper;            % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        options_idea2(1) = 400*minmax_problem.dim_d;          % number of function evaluations for IDEA
        options_idea2(4) = 5;                                 % number of agents
        options_idea2(31) = 1;                                % number of local restarts
        options_idea2(32) = 1;                                % communication heuristics (see com_red.m)
        options_idea2(22) = 0.1;                              % delta_c: crowding factor for global restart
        options_idea2(27) = 0.25;                             % tolconv: local convergence
        options_idea2(15) = -1;                               % verbousity
        options_idea2(28) = 0.2;                                % F
        options_idea2(29) = 0.8;                                % CR
        % options_idea2(30) = minmax_problem.n_obj;             % number of objectives (not used)
        options_idea2(33) = 0.2;                              % initial probability of running the subproblem (not used)

algo_outer.par.options = options_idea2;
algo_outer.par.initial_population = []; 


%% ALGORITHM INNER LOOP: IDEA2
algo_inner.optimise = @optimise_idea2_wrapper;            % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        options_idea2(1) = 150*minmax_problem.dim_u;          % number of function evaluations for IDEA
        options_idea2(4) = 5;                                 % number of agents
        options_idea2(31) = 1;                                % number of local restarts
        options_idea2(32) = 1;                                % communication heuristics (see com_red.m)
        options_idea2(22) = 0.1;                              % delta_c: crowding factor for global restart
        options_idea2(27) = 0.25;                             % tolconv: local convergence
        options_idea2(15) = -1;                               % verbousity
        options_idea2(28) = 0.2;                                % F
        options_idea2(29) = 0.8;                                % CR
        % options_idea2(30) = minmax_problem.n_obj;             % number of objectives (not used)
        options_idea2(33) = 0.2;                              % initial probability of running the subproblem (not used)

algo_inner.par.options = options_idea2;
algo_inner.par.initial_population = [];

return
