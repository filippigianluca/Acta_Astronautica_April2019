function [d,fval,exitflag,output] = optimise_relax_bisurr(problem_minmax, algo_outer, algo_inner, par_minmax)

global nfevalglobal
% Rename inputs
n_d = problem_minmax.dim_d;
n_u = problem_minmax.dim_u;
n_obj = problem_minmax.n_obj;
sign_inner = problem_minmax.sign_inner; % 1 for minmax, -1 for minmin

nfevalmax = par_minmax.maxnfeval;

% Build metaproblems
problem_max_u = build_metaproblem_macsminmax_inner(problem_minmax);
problem_min_d = build_metaproblem_mo_sa_outer(problem_minmax);

% Random initial guess for d
n_d0 = par_minmax.n_d0;
% d_0 = lhsu(zeros(1,n_d),ones(1,n_d),n_d0);
d_0 = lhsdesign(n_d0,n_d,'criterion','maximin');

convergence = []; % record to watch how it converges

% Initialise archive
record = nan(n_d0,n_d+n_obj*(n_u+1));

for i = 1:size(d_0,1)
    problem_max_u.par_objfun.d = d_0(i,:);                                                 % tell metaproblem what d to fix
    u0i = [];
    fi = [];
    for obj = 1:n_obj
        problem_max_u.par_objfun.objective = obj;                                                       % tell metaproblem what objective to optimise
        [ umax, f_aux, ~ , output_aux ] = algo_inner{obj}.optimise(problem_max_u,algo_inner{obj}.par);  % optimise
        u0i = [u0i, umax];
        fi = [fi, -sign_inner*f_aux];
    end
    record(i,:) = [d_0(i,:), u0i, fi];
end

% initialise the reference and the record_aux
record_aux = record;
ref = record(dominance(record(:,end-n_obj+1:end),0) == 0,:);
ref 


% Main loop
iter = 0;
stop = false;
while ~stop
    iter=iter+1;

    % figure(14)
    % plot3(record(:,1),record(:,2),record(:,end),'.')
    % drawnow

    % manage dataset
    points_total = size(record,1);
    if( iter == 1 || points_total<= problem_min_d.par_objfun.surrogate.max_set_size)
        dataset = record;                                   % nothing to manage
    else
        dataset = record_aux;                               % start by adding the guys from last iteration

        record_data = record;                               % copy the record
        dom = dominance(record_data(:,end-n_obj+1:end),0);
        [dom,idx_sort] = sort(dom);
        record_data = record_data(idx_sort,:);              % the record is now sorted by dominance and dom is corresponding
        available = ~ismember(record_data,dataset,'rows');  % these are the guys you can sample
        
        best_guys = record_data(available(1:problem_min_d.par_objfun.surrogate.set_size_minima),:); % in MO there should be a shrink here and more ingenuity
                                                                                                    % for the moment just take the best you have, IF available    
        dataset=[dataset;best_guys];
        available(1:problem_min_d.par_objfun.surrogate.set_size_minima) = 0;    % these guys are no longer available, were they or not

        % at this point is possible that you have exceeded the surrogate size if your paremeters were wrong. in that case you datasample the dataset
        datasize_1 = size(dataset,1);
        if (datasize_1 >= problem_min_d.par_objfun.surrogate.max_set_size)
            dataset = datasample(dataset,problem_min_d.par_objfun.surrogate.max_set_size,1,'Replace',false);
            warning('set_size_minima + last_iteration_additions >= max_set_size in surrogate building');
        elseif (sum(available)>0)
            % still points to sample, how many?
            remaining = problem_min_d.par_objfun.surrogate.max_set_size - datasize_1;
            if (remaining <= sum(available))
                dataset = [dataset ; datasample(record_data(available,:),remaining,1,'Replace',false)];
            else
                dataset = [dataset; record_data(available,:)];
            end
        end
    end

    % train surrogate
    problem_min_d.par_objfun.surrogate.model = problem_min_d.par_objfun.surrogate.training(dataset(:,1:n_d),dataset(:,n_d+1:end),problem_min_d.par_objfun.surrogate);   

    % "OUTER" LOOP: Compute dmin(i) = arg min {max f1(d,ue1), max f2(d,ue2)}
    [dmin, f_outer, ~ , output_aux] = algo_outer.optimise(problem_min_d,algo_outer.par);

    dmin(dmin < 0) = 0;
    dmin(dmin > 1) = 1;

    % Remove solutions returned more than once (if any)
    [~,idx] = unique(round(1e8*dmin),'rows');
    dmin=dmin(idx,:);
    f_outer = f_outer(idx,:);

    if (n_obj == 2)
        % figure(2)
        colors = 'ckbrmg';
        hold on
        plot(f_outer(:,1),f_outer(:,2),strcat(colors(mod(iter,length(colors))+1),'.'))
        drawnow
        % figure(1)    
    end
    d_outer = dmin.*repmat(problem_minmax.ub_d'-problem_minmax.lb_d',[size(dmin,1),1]) + repmat(problem_minmax.lb_d',[size(dmin,1),1])
    f_outer

    % "INNER" LOOP: Compute ue(i+1) = arg max f(dmin(i),u)
    record_aux = [];
    for i = 1:size(dmin,1)
        problem_max_u.par_objfun.d = dmin(i,:);
        ufmax = nan(1,(n_u+1)*n_obj);
        [uf0,~] = problem_min_d.par_objfun.surrogate.predictor(dmin(i,:), problem_min_d.par_objfun.surrogate.model);

        for obj = 1:n_obj
            problem_max_u.par_objfun.objective = obj;

            % Give a member of the population by psi-surr
            algo_inner{obj}.par.initial_population=uf0((obj-1)*n_u+1 : obj*n_u);

            % Maximize subproblem to find umax = arg max f(dmin(i),u)
            [ umax, f_inner , ~ , output_aux] = algo_inner{obj}.optimise(problem_max_u,algo_inner{obj}.par);
            % nfeval = nfeval + output_aux.nfeval;
            umax(umax < 0) = 0;
            umax(umax > 1) = 1;
            f_inner = -sign_inner*f_inner;
            
            % ALWAYS ARCHIVE
            ufmax((obj-1)*n_u+1 : obj*n_u) = umax;
            ufmax(n_u*n_obj+obj) = f_inner;
        end

        % archive for dmin
        record_aux = [record_aux ; dmin(i,:), ufmax];
    end
    
    if (n_obj == 2)
        % figure(2)
        colors = 'ckbrmg';
        hold on
        plot(record_aux(:,end-1),record_aux(:,end),strcat(colors(mod(iter,length(colors))+1),'o'))
        drawnow
        % figure(1)    
    end
    f_inner = record_aux(:,end-n_obj+1:end)

    % Archive dmin
    record = [record_aux ; record]; % concatenate in the beginning so that unique after will priorize

    %% I THINK THE CROSS CHECKS SHOULD BE HERE IN CASE YOU CONVERGED TO A GUY YOU ALREADY KNOW
    
    % Remove solutions archived more than once (if any)
    [~, idx] = unique(round(1e8*record(:,1:n_d)),'rows','stable');
    record = record(idx,:);

    % update the reference
    ref = record(dominance(record(:,end-n_obj+1:end),0) == 0 ,:);
    ref
    nfevalglobal

    convergence = [convergence; nfevalglobal.*ones(size(ref,1)), ref(:,end-n_obj+1:end)];
    stop = nfevalglobal >= nfevalmax;
end

fval = ref(:,end-n_obj+1:end);
frontsize = size(ref,1);

d = ref(:,1:n_d).*repmat(problem_minmax.ub_d'-problem_minmax.lb_d',[frontsize,1]) + repmat(problem_minmax.lb_d',[frontsize,1]);
u = cell(1,n_obj);
for obj = 1:n_obj
    map_info = problem_max_u.par_objfun.map_u_info{obj};
    for i = 1:frontsize
        u{obj}(i,:) = map_affine(ref(i,n_d+(obj-1)*n_u+1 : n_d+obj*n_u),map_info); %this can be easily vectorized
    end
end

output.u = u;
output.nfeval = nfevalglobal;
output.iter = iter;
output.convergence = convergence;
output.record = record;
exitflag = 0;

end