function [bel] = sstl2p(design)
bel = 0.0;
objective = 1;
n_iter = 7;
presample = 200;
opt_count = 0;

init = str2func(strcat('init_algo_so_sstl_corners'));
savefolder = strcat('RESULTS/SC/');

% global nfevalglobal;
% nfevalglobal = 0;

% global f_history;
% f_history=[];

%% initialise problem
[ problem_0 ] = init_problem_sstl2();
problem_0.split_coordinates = [];

dim_u = problem_0.dim_u;
lb_u = problem_0.lb_u{objective};
ub_u = problem_0.ub_u{objective};
problem_list_next = [problem_0];
% n_int_orig = cellfun('size',problem_0.lb_u{objective},2);

[ ~, ~, algo_inner ] = init(problem_0);
problem_max_u = build_metaproblem_macsminmax_inner(problem_0);
problem_max_u.par_objfun.objective = objective;
problem_max_u.par_objfun.d = (design-problem_0.lb_d')./(problem_0.ub_d'-problem_0.lb_d');

%% optimise
[ umax, fmax , ~ , output_aux] = algo_inner.optimise(problem_max_u,algo_inner.par);
opt_count = opt_count+1;
fmax = -fmax;
problem_0.provisional_fmax = fmax;

if(fmax<=0)
    bel = 1.0;
else
    [ umin, fmin , ~ , output_aux] = algo_inner.minimise(problem_max_u,algo_inner.par);
    fmin = -fmin;
    
    if (fmin<0) %otherwise nothing to do but exit and return bel = 0
        n_int = cellfun('size',problem_0.lb_u{objective},2);
        bpa = problem_0.bpa{objective};
        theta_list = [];
        
        % presample the theta structure for first analysis
        if (presample > 0)
            fe_pre = presample_fe(presample,n_int);
            for i = 1:presample
                problem_1 = problem_0;
                problem_1_bpa = 1;
                for j=1:dim_u
                    
                    problem_1.lb_u{objective}{j} = problem_0.lb_u{objective}{j}(1,fe_pre(i,j));
                    problem_1.ub_u{objective}{j} = problem_0.ub_u{objective}{j}(1,fe_pre(i,j));
                    problem_1_bpa = problem_1_bpa * problem_0.bpa{objective}{j}(1,fe_pre(i,j));
                end
                
                problem_max_u = build_metaproblem_macsminmax_inner(problem_1);
                problem_max_u.par_objfun.objective = objective;
                problem_max_u.par_objfun.d = (design-problem_0.lb_d')./(problem_0.ub_d'-problem_0.lb_d');
                
                [ umax, fmax , ~ , output_aux] = algo_inner.optimise(problem_max_u,algo_inner.par);
                fmax = -fmax;
                theta_list(end+1,:)=[pos2nfe(fe_pre(i,:),n_int), fe_pre(i,:) , fmax , problem_1_bpa];
            end
        end
        
        %             last_f_history = 0;
        problem_list_next = [problem_0];
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
                problem_max_u.par_objfun.d = (design-problem_minmax.lb_d')./(problem_minmax.ub_d'-problem_minmax.lb_d');
                
                
                %% optimise
                [ umax, fmax , ~ , output_aux] = algo_inner.optimise(problem_max_u,algo_inner.par);
                
                fmax= max(-fmax,problem_minmax.provisional_fmax);
                opt_count = opt_count+1;
                
                prob_maxima_list = [prob_maxima_list; fmax];
                prob_bpa = 1;
                for d = 1:dim_u
                    prob_bpa = prob_bpa * sum(problem_minmax.bpa{objective}{d});
                end
                prob_bpa_list = [prob_bpa_list; prob_bpa];
                
                u_true = map_affine(umax,problem_max_u.par_objfun.map_u_info{objective});
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
                    if isempty(find(theta_list(:,1)==n_fe_i))
                        problem_bpa = 1;
                        for j=1:dim_u
                            problem_bpa = problem_bpa * bpa{j}(position_fe(pos,j));
                        end
                        theta_list(end+1,:)=[pos2nfe(position_fe(pos,:),n_int), position_fe(pos,:) , fmax , problem_bpa];
                    end
                end
                
                
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
                if (isempty(maxima))
                    a=0;
                end
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
            
            % % save(strcat('SSTL_bel_iter_',num2str(iter)));
            % plot_belief(prob_maxima_list,prob_bpa_list);
            % title(num2str(iter));
        end
        [prob_maxima_list,idx] = sort(prob_maxima_list);
        prob_bpa_list = prob_bpa_list(idx);
        bel = sum(prob_bpa_list(prob_maxima_list<=0));
    end
    % plot_belief(prob_maxima_list,prob_bpa_list);
end
return


% %create results directory
% mkdir(savefolder);
% save(strcat(savefolder,'SC_',num2str(runid)));

% end