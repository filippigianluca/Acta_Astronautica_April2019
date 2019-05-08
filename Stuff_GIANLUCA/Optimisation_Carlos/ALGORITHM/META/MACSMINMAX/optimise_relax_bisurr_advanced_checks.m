function [d,fval,exitflag,output] = optimise_relax_bisurr_advanced_checks(problem_minmax, algo_outer, algo_inner, par_minmax)

global nfevalglobal
tol = 1e8; 
% Rename inputs
n_d = problem_minmax.dim_d;
n_u = problem_minmax.dim_u;
n_obj = problem_minmax.n_obj;
sign_inner = problem_minmax.sign_inner; % 1 for minmax, -1 for minmin

nfevalmax = par_minmax.maxnfeval;

% Cross-checks
cc_opt = par_minmax.cc_opt;
cc_opt.tol = tol;

% Build metaproblems
problem_max_u = build_metaproblem_macsminmax_inner(problem_minmax);
problem_min_d = build_metaproblem_mo_sa_outer(problem_minmax);

% Random initial guess for d
n_d0 = par_minmax.n_d0;
% d_0 = lhsu(zeros(1,n_d),ones(1,n_d),n_d0);
d_0 = lhsdesign(n_d0,n_d,'criterion','maximin'); % there shouldn't be repeated d's

convergence = []; % record to watch how it converges

% Initialise archive with maximisations
record = nan(n_d0,n_d+n_obj*(n_u+1));
d_info = struct;
u_record = cell(1,n_obj);
for i = 1:size(d_0,1)
    problem_max_u.par_objfun.d = d_0(i,:);                                                 % tell metaproblem what d to fix
    u0i = [];
    fi = [];
    
    d_info(i).d = d_0(i,:);
    d_info(i).umax = zeros(1,n_obj);
    d_info(i).fmax = nan(1,n_obj);
    d_info(i).u_checked = cell(1,n_obj);
    d_info(i).u_ls = cell(1,n_obj);

    for obj = 1:n_obj
        problem_max_u.par_objfun.objective = obj;                                                       % tell metaproblem what objective to optimise
        [ umax, f_aux, ~ , output_aux ] = algo_inner{obj}.optimise(problem_max_u,algo_inner{obj}.par);  % optimise
        u0i = [u0i, umax];
        fi = [fi, -sign_inner*f_aux];
        d_info(i).fmax(obj) = -sign_inner*f_aux;

        if isempty(u_record{obj})
            loc = 0;
        else
            
            [~, loc] = ismember(round(tol*umax),round(tol*u_record{obj}),'rows');
        end
        if loc > 0
            d_info(i).u_checked{obj}(end+1) = loc;
            d_info(i).umax(obj) = loc;
        else
            u_record{obj}(end+1,:) = umax;
            l = size(u_record{obj},1);
            d_info(i).u_checked{obj}(end+1) = l;
            d_info(i).umax(obj) = l;
        end        
    end
    record(i,:) = [d_0(i,:), u0i, fi];
end

% Do the first cross-checks
[record, d_info, u_record, nfeval_cc] = cross_check(problem_max_u, record, d_info,u_record,cc_opt);


% initialise the reference and the record_aux
record_aux = record;
size_record_aux = size(record_aux,1);
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
        dataset = record(end-size_record_aux+1:end,:);      % start by adding the guys from last iteration

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
    [~,idx] = unique(round(tol*dmin),'rows');
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
    dmin
    f_outer

    % "INNER" LOOP: Compute ue(i+1) = arg max f(dmin(i),u)
    record_aux = [];
    d_info_aux = [];
    for i = 1:size(dmin,1)
        % even if you already know that d, you redo the maximisation bc you converged there
        di = dmin(i,:);
        problem_max_u.par_objfun.d = di;
        ufmax = nan(1,(n_u+1)*n_obj);
        [uf0,~] = problem_min_d.par_objfun.surrogate.predictor(di, problem_min_d.par_objfun.surrogate.model);

        for obj = 1:n_obj
            problem_max_u.par_objfun.objective = obj;

            % Give a member of the population by psi-surr
            u0 = uf0((obj-1)*n_u+1 : obj*n_u);
            u0(u0 < 0) = 0;
            u0(u0 > 1) = 1;
            algo_inner{obj}.par.initial_population= u0;

            % Maximize subproblem to find umax = arg max f(dmin(i),u)
            [ umax, f_inner , ~ , output_aux] = algo_inner{obj}.optimise(problem_max_u,algo_inner{obj}.par);
            % nfeval = nfeval + output_aux.nfeval;
            umax(umax < 0) = 0;
            umax(umax > 1) = 1;
            f_inner = -sign_inner*f_inner;
            
            % remember u and fmax
            ufmax((obj-1)*n_u+1 : obj*n_u) = umax;
            ufmax(n_u*n_obj+obj) = f_inner;
        end

        % archive for dmin
        [~, loc] = ismember(round(tol*di),round(tol*record(:,1:n_d)),'rows');
        if loc>0 % you have converged to a d you already knew and rerun the maximisation
            % NOTE: in this case we are treating the d's as equal, since they are equal to the tolerance set.
            %% IF they are actually not regarding objective function, it's a matter of tolerance and its user's fault
            % so let's check if the new maximisations did something
            for obj = 1:n_obj
                fiobj = ufmax(n_u*n_obj+obj);
                uiobj = ufmax((obj-1)*n_u+1 : obj*n_u);
                [~, loc2] = ismember(round(tol*uiobj),round(tol*u_record{obj}),'rows');
                if loc2 > 0 % you already have this u in the archive. no need to add it
                    % but maybe you had never checked this d and this u together
                    if ~ismember(loc2,d_info(loc).u_checked{obj})
                        d_info(loc).u_checked{obj}(end+1) = loc2;
                    end
                    % and maybe it maximises
                    if sign_inner*fiobj > sign_inner*record(loc,n_d+n_u*n_obj+obj)
                        record(loc , n_d+(obj-1)*n_u+1 :n_d+ obj*n_u) = uiobj;
                        record(loc , n_d+n_u*n_obj+obj) = fiobj;
                        d_info(loc).umax(obj) = loc2;
                        d_info(loc).fmax(obj) = fiobj;
                    end

                else % there is a new u associated
                    u_record{obj}(end+1,:) = uiobj;     % archive it
                    l = size(u_record{obj},1);       
                    d_info(loc).u_checked{obj}(end+1) = l;   % remember it has already been cross-checked
                    % now does it maximise?
                    if sign_inner*fiobj > sign_inner*record(loc,n_d+n_u*n_obj+obj)
                        record(loc , n_d+(obj-1)*n_u+1 :n_d+ obj*n_u) = uiobj;
                        record(loc , n_d+n_u*n_obj+obj) = fiobj;
                        d_info(loc).umax(obj) = l;
                        d_info(loc).fmax(obj) = fiobj;
                    end                   
                end
            end

        else % this is a new d that has never been cross-checked

            record_aux = [record_aux ; di, ufmax]; % you add it to record_aux with its current u and fmax, since it was never cross-checked
            
            % you initialise the infos
            di_info = struct;
            di_info.d = di;
            di_info.umax = zeros(1,n_obj);
            di_info.fmax = nan(1,n_obj);
            di_info.u_checked = cell(1,n_obj);
            di_info.u_ls = cell(1,n_obj);

            for obj = 1:n_obj
                di_info.fmax(obj) = ufmax(n_u*n_obj+obj);
                uiobj = ufmax((obj-1)*n_u+1 : obj*n_u);
                [~, loc2] = ismember(round(tol*uiobj),round(tol*u_record{obj}),'rows');
                if loc2 > 0 % you already have this u in the archive
                    di_info.u_checked{obj}(end+1) = loc2;
                    di_info.umax(obj) = loc2;
                else % the new d has a new u associated
                    u_record{obj}(end+1,:) = umax;
                    l = size(u_record{obj},1);
                    di_info.u_checked{obj}(end+1) = l;
                    di_info.umax(obj) = l;
                end
            end

            d_info_aux = [d_info_aux, di_info]; % and you also remember his infos
        end
    end
    
    if (n_obj == 2)
        % figure(2)
        colors = 'ckbrmg';
        hold on
        plot(record_aux(:,end-1),record_aux(:,end),strcat(colors(mod(iter,length(colors))+1),'o'))
        drawnow
        % figure(1)    
    end
    
    if ~isempty(record_aux)
        f_inner = record_aux(:,end-n_obj+1:end)
    end
    % Archive dmin and infos
    size_record_aux = size(record_aux,1); % will come in handy for next iteration
    record = [record; record_aux];
    d_info = [d_info, d_info_aux];

    %% I THINK THE CROSS CHECKS SHOULD BE HERE IN CASE YOU CONVERGED TO A GUY YOU ALREADY KNOW
    [record, d_info, u_record, nfeval_cc] = cross_check(problem_max_u, record, d_info,u_record,cc_opt);
    
    % % Remove solutions archived more than once (if any) --> THIS SHOULDN'T BE NECESSARY BC WE TOOK CARE IN ARCHIVING
    % [~, idx] = unique(round(tol*record(:,1:n_d)),'rows','stable');
    % record = record(idx,:);

    % update the reference
    ref = record(dominance(record(:,end-n_obj+1:end),0) == 0 ,:);
    ref
    nfevalglobal
    if (n_obj == 2)
        title(num2str(nfevalglobal))
        drawnow
    end
    % figure(1)
    % plot(record(:,1),record(:,3),'.')
    % drawnow
    % figure(2)
    % plot(record(:,1),record(:,2),'.')
    % drawnow
    
    
    convergence = [convergence; nfevalglobal.*ones(size(ref,1),1), ref(:,end-n_obj+1:end)];
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
output.u_record = u_record;
output.nfeval = nfevalglobal;
output.iter = iter;
output.convergence = convergence;
output.record = record;
exitflag = 0;

end