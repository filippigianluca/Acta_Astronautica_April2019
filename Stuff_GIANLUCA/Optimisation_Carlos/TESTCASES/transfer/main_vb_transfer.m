%% REMEMBER TO ACTIVATE ARCHIVING IN MASK MACSMINMAX INNER

% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_transfer_vb'));
savefolder = strcat('RESULTS/TRANSFER_BEL/');

dminmax = [ 62363.1505101023   17.8344027493367 ];
fminmax = [ 22.8344027493367    0.877322119681857];
reference_values = [fminmax(1,1)*1.1   10^fminmax(1,2)*1.1]; %for computing the score AD HOC for this problem
scale_factors = [500 250];

% objective = 1;
n_iter = 7;
presample = 0;

    global nfevalglobal;
    nfevalglobal = 0;
    
    global f_history;
    % global theta_history;
    f_history=[];
    % theta_history =[];

    opt_count = 0;
    
    %% initialise problem
    [ problem_0 ] = init_problem_transfer_vb();
    problem_0.split_coordinates = [];
    problem_0.provisional_fmax = nan(1,problem_0.n_obj);

    problem_list_next = [problem_0];
    % n_int_orig = cellfun('size',problem_0.lb_u{objective},2);
    
    %% build structure for approx belief curve
    % NOTE that I AM TAKING ONLY ONE U-SSTRUCTURE FOR ALL OBJECTIVES
    lb_u = problem_0.lb_u{1,1}; 
    ub_u = problem_0.ub_u{1,1};
    bpa = problem_0.bpa{1,1};
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

    %% FINISH PRESAMPLE!
    % % for cases like sstl where you maximise with just 1 feval, you need to presample
    % if (presample > 0)
    %     u_sc = lhsgen(presample, dim_u);
    %     for i = 1:presample
    %         problem_max_u = build_metaproblem_macsminmax_inner(problem_minmax);
    %         u_true = map_affine(u_sc(i,:),problem_max_u.par_objfun.map_u_info{objective});


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
            
            problem_max_u.par_objfun.d = (dminmax-problem_minmax.lb_d')./(problem_minmax.ub_d'-problem_minmax.lb_d');


            %% optimise
            fmax = [];
            for obj = 1:problem_minmax.n_obj
                problem_max_u.par_objfun.objective = obj;
                [ umax_obj, fmax_obj , ~ , output_aux] = algo_inner{obj}.optimise(problem_max_u,algo_inner{obj}.par);

                fmax_obj= max(-fmax_obj,problem_minmax.provisional_fmax(1,obj));
                opt_count = opt_count+1;
                fmax(1,obj) = fmax_obj;
            end
            % [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

            prob_maxima_list = [prob_maxima_list; fmax];
            prob_bpa = 1;
            for d = 1:dim_u
                prob_bpa = prob_bpa * sum(problem_minmax.bpa{1,1}{d});
            end
            prob_bpa_list = [prob_bpa_list; prob_bpa];

            %% update approx belief ONLY FOR deltaV
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

            last_f_history = 0;
            f_history = [];


            theta_list_aux = theta_list;
            for num_split = 1:size(problem_minmax.split_coordinates,1)
                theta_list_aux = theta_list_aux(theta_list_aux(:,1+problem_minmax.split_coordinates(num_split,1)) ==  problem_minmax.split_coordinates(num_split,2) ,:);
            end
            next_dim_split = 0;
            score_max_tot = -inf;
            for d=1:dim_u
                if isempty(problem_minmax.split_coordinates) || ~any(problem_minmax.split_coordinates(:,1)==d)
                    
                    dvlog_maxima = [];
                    tof_maxima = [];
                    bpas = [];
                    for i_int = 1:n_int(d)
                        % find dv maxima with the history and tof maxima with an evaluation
                        theta_list_aux_aux = theta_list_aux(theta_list_aux(:,1+d)==i_int,:);
                        if (~isempty(theta_list_aux_aux))
                            dvlog_maxima = [dvlog_maxima; max(theta_list_aux_aux(:,end-1))];
                            bpas = [bpas; bpa{d}(i_int)];

                            problem_aux_tof = problem_minmax;
                            problem_aux_tof.lb_u{2,1}{d,1} = lb_u{d};
                            problem_aux_tof.ub_u{2,1}{d,1} = ub_u{d};

                            problem_aux_tof_max_u = build_metaproblem_macsminmax_inner(problem_minmax);
            
                            problem_aux_tof_max_u.par_objfun.d = (dminmax-problem_minmax.lb_d')./(problem_minmax.ub_d'-problem_minmax.lb_d');
                            problem_aux_tof_max_u.par_objfun.objective = 1;
                            [ umax_aux_tof, fmax_aux_tof , ~ , output_aux] = algo_inner{1}.optimise(problem_aux_tof_max_u,algo_inner{1}.par);
                            tof_maxima = [tof_maxima; -fmax_aux_tof];

                        end
                    end

                    %% COMPUTE SCORE
                    dv_maxima = 10.^dvlog_maxima;
                    score = 0;
                    for k=1:length(bpas)
                        score = score + abs((fmax(1)-tof_maxima(k))*(10^fmax(2)-dv_maxima(k)));
                        score = score + abs((reference_values(1)-fmax(1))*(10^fmax(2)-dv_maxima(k)));
                        score = score + abs((reference_values(2)-10^fmax(2))*(fmax(1)-tof_maxima(k)));
                        score = score * bpas(k);
                    end

                    if score>score_max_tot
                        next_dim_split = d;
                        score_max_tot = score;
                    end
                    % d
                    % [tof_maxima, dv_maxima]
                    % score
                end
            end
            iter
            next_dim_split
            % problem_list_next_aux = split_problem(problem_minmax,next_dim_split);
            % umax_true = map_affine(umax,problem_max_u.par_objfun.map_u_info{objective});
            problem_list_next_aux = [];
            for i_prob = 1:n_int(next_dim_split)
                problem_aux = problem_minmax;
                for obj=1:problem_aux.n_obj
                    problem_aux.lb_u{obj}{next_dim_split} = problem_aux.lb_u{obj}{next_dim_split}(1,i_prob);
                    problem_aux.ub_u{obj}{next_dim_split} = problem_aux.ub_u{obj}{next_dim_split}(1,i_prob);
                    problem_aux.bpa{obj}{next_dim_split} = problem_aux.bpa{obj}{next_dim_split}(1,i_prob);
                end
                problem_aux.split_coordinates = [problem_aux.split_coordinates; next_dim_split, i_prob];
                % if problem_aux.lb_u{objective}{next_dim_split} <= umax_true(next_dim_split) && umax_true(next_dim_split) <= problem_aux.ub_u{objective}{next_dim_split}
                %     problem_aux.provisional_fmax = fmax;
                % else
                %     problem_aux.provisional_fmax = nan;
                % end
                problem_list_next_aux = [problem_list_next_aux; problem_aux];
            end
            
            problem_list_next = [problem_list_next; problem_list_next_aux];
        
        end

        prob_maxima_list     
        xx = 17:0.05:reference_values(1);
        yy = 5.5 : 0.005: reference_values(2);
        [XX,YY] = meshgrid(xx,yy);
        zz = zeros(size(XX));
        for ix=1:length(xx)
            for jy = 1:length(yy)
                for m = 1:size(prob_maxima_list,1)
                    if xx(ix)>=prob_maxima_list(m,1) && yy(jy)>=10^prob_maxima_list(m,2)

                        zz(jy,ix)=zz(jy,ix)+prob_bpa_list(m);
                    end
                end
            end
        end
        figure(1)
        mesh(XX,YY,zz)
        drawnow
        xlabel('tof')
        ylabel('dv')
        zlabel('Bel')
        title(strcat('iter ',num2str(iter),' : ', num2str(size(prob_maxima_list,1)), ' FE'))
        % %create results directory
        mkdir(savefolder);
        save(strcat(savefolder,'transfer_bel_fe_',num2str(size(prob_maxima_list,1))));
    end


% end