%% REMEMBER TO ACTIVATE ARCHIVING IN MASK MACSMINMAX INNER

% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_so'));
savefolder = strcat('RESULTS/SC/');

% dminmax = [10  90  7   0.380763886588873   0.166313359938401   0.278354987823124   145 3   1   0.750742712651230];
dminmax = [10.0454441443253    89.7707372551649    7.00017189812877    0.465883837532156   0.330767275050868   0.296528323438445   144.696780602974    3.00065175047948    1.00274021888474    0.830477346032066];
objective = 1;
n_iter = 16;
% dminmax = [10.0     90.0   7   0.25  1.0  0.28    145     3.0     1.0     0.75 ];
%%
% for runid= 1     % 20 runs
% 
    % disp(strcat('SC_',num2str(runid)))
    global nfevalglobal;
    nfevalglobal = 0;
    
    global f_history;
    % global theta_history;
    f_history=[];
    % theta_history =[];

    opt_count = 0;
    
    %% initialise problem
    [ problem_0 ] = init_problem_sc();
    problem_0.split_coordinates = [];
    problem_0.provisional_fmax = nan;

    problem_list_next = [problem_0];
    % n_int_orig = cellfun('size',problem_0.lb_u{objective},2);
    
    %% build structure for approx belief curve
    lb_u = problem_0.lb_u{objective}; 
    ub_u = problem_0.ub_u{objective};
    bpa = problem_0.bpa{objective};
    dim_u = problem_0.dim_u;
    n_int = cellfun('size',lb_u,2);
    n_fe_tot = prod(n_int);
    theta_list = [(1:n_fe_tot)',nan(n_fe_tot,dim_u+1), ones(n_fe_tot,1)];
    for i=1:n_fe_tot
        pos = nfe2pos(i,n_int);
        theta_list(i,2:dim_u+1) = pos;
        for j = 1:length(pos)
            theta_list(i,end) = theta_list(i,end) * bpa{j}(pos(j));
        end
    end
    last_f_history = 0;


    for iter = 1:n_iter
        problem_list = problem_list_next;
        problem_list_next = [];

        prob_maxima_list = [];
        prob_bpa_list = [];
        for problem_id = 1:length(problem_list)
            problem_minmax = problem_list(problem_id);
            %% initialise algorithm
            [ ~, ~, algo_inner ] = init(problem_minmax);
            problem_max_u = build_metaproblem_macsminmax_inner(problem_minmax);
            problem_max_u.par_objfun.objective = objective;
            problem_max_u.par_objfun.d = (dminmax-problem_minmax.lb_d')./(problem_minmax.ub_d'-problem_minmax.lb_d');


            %% optimise
            [ umax, fmax , ~ , output_aux] = algo_inner.optimise(problem_max_u,algo_inner.par);

            fmax= max(-fmax,problem_minmax.provisional_fmax);
            opt_count = opt_count+1;
            % [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

            prob_maxima_list = [prob_maxima_list; fmax];
            prob_bpa = 1;
            for d = 1:dim_u
                prob_bpa = prob_bpa * sum(problem_minmax.bpa{objective}{d});
            end
            prob_bpa_list = [prob_bpa_list; prob_bpa];

            %% update approx belief
            global f_history;
            for i= last_f_history+1:size(f_history,1)
                % u = f_history(i,1:end-1);
                % u_true = map_affine(u,problem_max_u.par_objfun.map_u_info{objective});
                u_true = f_history(i,1:end-1);
                position_fe = zeros(1,dim_u);
                for d = 1:dim_u
                    for i_int = 1:length(lb_u{d})
                        if lb_u{d}(1,i_int) <= u_true(d) && u_true(d) <= ub_u{d}(1,i_int)
                            if position_fe(1,d) == 0
                                position_fe(:,d) = i_int;
                            else
                                position_fe_aux = position_fe;
                                position_fe_aux(:,d) = i_int;
                                position_fe = [position_fe; position_fe_aux];
                            end
                        end
                    end
                end

                for pos = 1:size(position_fe,1)
                    n_fe_i = pos2nfe(position_fe(pos,:), n_int);
                    if isnan(theta_list(n_fe_i,end-1)) || f_history(i,end) > theta_list(n_fe_i,end-1)
                        theta_list(n_fe_i,end-1) = f_history(i,end);
                    end
                end
            end

            % last_f_history = i;

            % last_f_history = 0;
            f_history = [];


            theta_list_aux = theta_list;
            for num_split = 1:size(problem_minmax.split_coordinates,1)
                theta_list_aux = theta_list_aux(theta_list_aux(:,1+problem_minmax.split_coordinates(num_split,1)) ==  problem_minmax.split_coordinates(num_split,2) ,:);
            end
            next_dim_split = 0;
            min_tot = inf;
            for d=1:dim_u
                if isempty(problem_minmax.split_coordinates) || ~any(problem_minmax.split_coordinates(:,1)==d)
                    maxima = [];
                    bpas = [];
                    for i_int = 1:n_int(d)
                        theta_list_aux_aux = theta_list_aux(theta_list_aux(:,1+d)==i_int,:);
                        if (~isempty(theta_list_aux_aux))
                            maxima = [maxima; max(theta_list_aux_aux(:,end-1))];
                            bpas = [bpas; bpa{d}(i_int)];
                        end
                    end

                    if (min(maxima)) < min_tot
                        min_tot = min(maxima);
                        next_dim_split = d; 
                    end
                      % if (~isempty(maxima))
                      %     plot_belief(maxima,bpas);
                      %     title(num2str(d))
                      % end
                end
            end
            
            % problem_list_next_aux = split_problem(problem_minmax,next_dim_split);
            umax_true = map_affine(umax,problem_max_u.par_objfun.map_u_info{objective});
            problem_list_next_aux = [];
            for i_prob = 1:n_int(next_dim_split)
                problem_aux = problem_minmax;
                problem_aux.lb_u{objective}{next_dim_split} = problem_aux.lb_u{objective}{next_dim_split}(1,i_prob);
                problem_aux.ub_u{objective}{next_dim_split} = problem_aux.ub_u{objective}{next_dim_split}(1,i_prob);
                problem_aux.bpa{objective}{next_dim_split} = problem_aux.bpa{objective}{next_dim_split}(1,i_prob);
                problem_aux.split_coordinates = [problem_aux.split_coordinates; next_dim_split, i_prob];
                if problem_aux.lb_u{objective}{next_dim_split} <= umax_true(next_dim_split) && umax_true(next_dim_split) <= problem_aux.ub_u{objective}{next_dim_split}
                    problem_aux.provisional_fmax = fmax;
                else
                    problem_aux.provisional_fmax = nan;
                end
                problem_list_next_aux = [problem_list_next_aux; problem_aux];
            end
            
            problem_list_next = [problem_list_next; problem_list_next_aux];
        
        end

        save(strcat('sc_bel_iter_',num2str(iter)));
        plot_belief(prob_maxima_list,prob_bpa_list);
        title(num2str(iter));
    end

    
    % %create results directory
    % mkdir(savefolder);
    % save(strcat(savefolder,'SC_',num2str(runid)));

% end