function [d,fval,exitflag,output] = optimise_ausa_trisurr(problem_minmax, algo_outer, algo_inner, par_minmax)

global nfevalglobal
tol = 1e8; 
% Rename inputs
n_d = problem_minmax.dim_d;
n_u = problem_minmax.dim_u;
n_obj = problem_minmax.n_obj;
sign_inner = problem_minmax.sign_inner; % 1 for minmax, -1 for minmin

nfevalmax = par_minmax.maxnfeval;
max_iter = inf;
if isfield(par_minmax,'max_iter')
    max_iter=par_minmax.max_iter;
end

verbosity = true;
if isfield(par_minmax,'verbosity')
    verbosity=par_minmax.verbosity;
end

% Cross-checks
cc_opt = par_minmax.cc_opt;
cc_opt.tol = tol;

% Build metaproblems and such
problem_max_u = build_metaproblem_macsminmax_inner(problem_minmax);
problem_eval = problem_max_u;

% for alternation
alternate_outer = par_minmax.alternate_outer;
meth_outer = [];
i_meth_outer = 0;
iter_alternation = 0;
if alternate_outer
    if algo_outer.par_sa.use
        meth_outer = [meth_outer, 1];
    end
    if algo_outer.par_sp.use
        meth_outer = [meth_outer, 2];
    end
    if algo_outer.par_au.use
        meth_outer = [meth_outer, 3];
    end
    if isempty(meth_outer)
        error('alternating between no methods')
    end
    i_meth_outer = 1;
end

subproblem_minmax = []; subalgo_minmax = []; subalgo_outer = []; subalgo_inner = [];
if algo_outer.par_sp.use
    subproblem_minmax = algo_outer.par_sp.init_subproblem(problem_minmax);
    [ subalgo_minmax, subalgo_outer, subalgo_inner ] = algo_outer.par_sp.init_subalgo(subproblem_minmax);
end

problem_min_d_au = [];
if algo_outer.par_au.use
    problem_min_d_au = build_metaproblem_macsminmax_outer(problem_max_u,algo_outer.par_au.local_search_flag);
end

psi_surr = false(1,n_obj);
for obj=1:n_obj
    psi_surr(obj)=algo_inner{obj}.par.psi_surr;
end

problem_min_d_sa = [];
if algo_outer.par_sa.use ||  any(psi_surr)
    problem_min_d_sa = build_metaproblem_mo_sa_outer(problem_minmax);
end



% Random initial guess for d
n_d0 = par_minmax.n_d0;
% d_0 = lhsu(zeros(1,n_d),ones(1,n_d),n_d0);
d_0 = lhsdesign(n_d0,n_d,'criterion','maximin'); % there shouldn't be repeated d's

% Random sample to maintain for S(d,u) -> otherwise the surrogate could degenerate if all the maximisations converge good to a single u for instance
n_du0 = 0;
du_0 = [];
fdu_0 = [];
duf_0 = [];
if (algo_outer.par_sp.use && n_du0 > 0)
    n_du0 = max(0,algo_outer.par_sp.n_du0);
    du_0 = lhsdesign(n_du0,n_d+n_u,'criterion','maximin');
    fdu_0 = nan(n_du0,n_obj);
    for i=1:n_du0
        problem_eval.par_objfun.d = du_0(i,1:n_d);
        for obj = 1:n_obj
            problem_eval.par_objfun.objective = obj;
            [fdu_0_obj] = problem_eval.objfun(du_0(i,n_d+1:n_d+n_u), problem_eval.par_objfun);
            fdu_0(i,obj) = -sign_inner*fdu_0_obj;
        end
    end

    duf_0 = [du_0(:,1:n_d), repmat(du_0(:,n_d+1:n_d+n_u),[1,n_obj]), fdu_0];

end

convergence = struct; % record to watch how it converges

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
if verbosity
    current_ref = ref;
    current_ref(:,1:n_d) = ref(:,1:n_d).*repmat(problem_minmax.ub_d'-problem_minmax.lb_d',[size(ref,1),1]) + repmat(problem_minmax.lb_d',[size(ref,1),1]);
    for obj = 1:n_obj
        map_info = problem_max_u.par_objfun.map_u_info{obj};
        for i = 1:size(ref,1)
            current_ref(i,n_d+(obj-1)*n_u+1:n_d+obj*n_u) = map_affine(ref(i,n_d+(obj-1)*n_u+1 : n_d+obj*n_u),map_info); %this can be easily vectorized
        end
    end
    current_ref
end

last_ref = ref; % for alternation

convergence(1).iter = 0;
convergence(1).nfeval = nfevalglobal;
convergence(1).d = ref(:,1:n_d);
convergence(1).u = ref(:,n_d+1:end-n_obj);
convergence(1).f = ref(:,end-n_obj+1:end);

% Main loop
iter = 0;
stop = false;
while ~stop
    iter=iter+1;
    
    if verbosity && par_minmax.alternate_outer
        i_meth_outer
    end
    % figure(14)
    % plot3(record(:,1),record(:,2),record(:,end),'.')
    % drawnow

    % manage dataset for s(d,u)
    if (algo_outer.par_sp.use &&(~alternate_outer || meth_outer(i_meth_outer) == 2))
        n_u_ref_div = max(0,algo_outer.par_sp.n_u_ref_diversity);
        points_total = size(record,1);
        if( iter == 1 || points_total+n_du0+n_u_ref_div<= subproblem_minmax.max_set_size)
            dataset_sdu = record;                          % nothing to manage
        else
            dataset_sdu = record(end-size_record_aux+1:end,:);      % start by adding the guys from last iteration

            record_data = record;                               % copy the record
            dom = dominance(record_data(:,end-n_obj+1:end),0);
            [dom,idx_sort] = sort(dom);
            record_data = record_data(idx_sort,:);              % the record is now sorted by dominance and dom is corresponding
            available = ~ismember(record_data,dataset_sdu,'rows');  % these are the guys you can sample
            
            best_guys = record_data(available(1:subproblem_minmax.set_size_minima),:); % in MO there should be a shrink here and more ingenuity
                                                                                                        % for the moment just take the best you have, IF available    
            dataset_sdu=[dataset_sdu;best_guys];
            available(1:subproblem_minmax.set_size_minima) = 0;    % these guys are no longer available, were they or not

            % at this point is possible that you have exceeded the surrogate size if your paremeters were wrong. in that case you datasample the dataset
            datasize_1 = size(dataset_sdu,1);
            if (datasize_1 >= subproblem_minmax.max_set_size - n_du0-n_u_ref_div)
                dataset_sdu = datasample(dataset_sdu,subproblem_minmax.max_set_size - n_du0-n_u_ref_div,1,'Replace',false);
                warning('set_size_minima + last_iteration_additions >= max_set_size in surrogate building');
            elseif (sum(available)>0)
                % still points to sample, how many?
                remaining = subproblem_minmax.max_set_size - n_du0-n_u_ref_div - datasize_1;
                if (remaining <= sum(available))
                    dataset_sdu = [dataset_sdu ; datasample(record_data(available,:),remaining,1,'Replace',false)];
                else
                    dataset_sdu = [dataset_sdu; record_data(available,:)];
                end
            end
        end

        if n_du0 > 0
            dataset_sdu = [dataset_sdu; duf_0]; % append the initial sample
        end

        if n_u_ref_div > 0
            du_ref_div = lhsdesign(n_u_ref_div,n_d+n_u,'criterion','maximin');
            d_span = repmat(max(0,algo_outer.par_sp.size_d_box_ref_diversity), [1, n_d]);
            ref_x_u_ref_div = nan(n_u_ref_div,n_d+(n_u+1)*n_obj);

            for i = 1:n_u_ref_div
                % pick a i_ref at random from ref+record_aux
                if isempty(record_aux)
                    d_ref = datasample(ref(:,1:n_d),1,1,'Replace',true);
                else
                    
                    d_ref = datasample([ref(:,1:n_d);record_aux(:,1:n_d)],1,1,'Replace',true);
                end
                % pick a point in the neighborhood
                d_ref_div = d_ref - d_span + 2.0*(d_span.*du_ref_div(i,1:n_d));
                d_ref_div(d_ref_div < 0) = 0;
                d_ref_div(d_ref_div > 1) = 1;
                ref_x_u_ref_div(i,1:n_d+n_obj*n_u) = [d_ref_div, repmat(du_ref_div(i,n_d+1:n_d+n_u), [1,n_obj])];

                problem_eval.par_objfun.d = d_ref_div;
                for obj = 1:n_obj
                    problem_eval.par_objfun.objective = obj;
                    [f_div_i_obj] = problem_eval.objfun(ref_x_u_ref_div(i,n_d+(obj-1)*n_u+1:n_d+obj*n_u), problem_eval.par_objfun);
                    ref_x_u_ref_div(i,n_d+n_obj*n_u+obj) = -sign_inner*f_div_i_obj;
                end
            end
            dataset_sdu = [dataset_sdu; ref_x_u_ref_div]; % append the ref x u_ref_div
        end

        % train surrogate
        for obj=1:n_obj
            subproblem_minmax.par_objfun{obj}.surrogate.model = subproblem_minmax.par_objfun{obj}.surrogate.training([dataset_sdu(:,1:n_d),dataset_sdu(:,n_d+(obj-1)*n_u+1:n_d+obj*n_u)] ,dataset_sdu(:,n_d+n_u*n_obj+obj),subproblem_minmax.par_objfun{obj}.surrogate);
        end
    end

    % manage dataset for Psi and Phi
    if ((any(psi_surr) || algo_outer.par_sa.use) && (~alternate_outer || meth_outer(i_meth_outer) == 1 || meth_outer(i_meth_outer) == 3))
        
        points_total = size(record,1);
        if( iter == 1 || points_total<= problem_min_d_sa.par_objfun.surrogate.max_set_size)
            dataset = record;                                   % nothing to manage
        else
            dataset = [record(1:n_d0,:);...
                       record(end-size_record_aux+1:end,:)];    % start by adding the guys from last iteration and the guys in n_d0

            record_data = record;                               % copy the record
            dom = dominance(record_data(:,end-n_obj+1:end),0);
            [dom,idx_sort] = sort(dom);
            record_data = record_data(idx_sort,:);              % the record is now sorted by dominance and dom is corresponding
            available = ~ismember(record_data,dataset,'rows');  % these are the guys you can sample
            
            best_guys = record_data(available(1:problem_min_d_sa.par_objfun.surrogate.set_size_minima),:); % in MO there should be a shrink here and more ingenuity
                                                                                                        % for the moment just take the best you have, IF available    
            dataset=[dataset;best_guys];
            available(1:problem_min_d_sa.par_objfun.surrogate.set_size_minima) = 0;    % these guys are no longer available, were they or not

            % at this point is possible that you have exceeded the surrogate size if your paremeters were wrong. in that case you datasample the dataset
            datasize_1 = size(dataset,1);
            if (datasize_1 >= problem_min_d_sa.par_objfun.surrogate.max_set_size)
                dataset = datasample(dataset,problem_min_d_sa.par_objfun.surrogate.max_set_size,1,'Replace',false);
                warning('set_size_minima + last_iteration_additions >= max_set_size in surrogate building');
            elseif (sum(available)>0)
                % still points to sample, how many?
                remaining = problem_min_d_sa.par_objfun.surrogate.max_set_size - datasize_1;
                if (remaining <= sum(available))
                    dataset = [dataset ; datasample(record_data(available,:),remaining,1,'Replace',false)];
                else
                    dataset = [dataset; record_data(available,:)];
                end
            end
        end

        % train surrogate
        problem_min_d_sa.par_objfun.surrogate.model = problem_min_d_sa.par_objfun.surrogate.training(dataset(:,1:n_d),dataset(:,n_d+1:end),problem_min_d_sa.par_objfun.surrogate);   

    end

    % "OUTER" LOOP with the min-max subproblem
    dmin_sp = []; f_outer_sp = []; u_outer_sp = cell(1,n_obj);
    if(algo_outer.par_sp.use && (~alternate_outer || meth_outer(i_meth_outer) == 2))
        [dmin_sp, f_outer_sp, ~ , output_aux] = subalgo_minmax.optimise(subproblem_minmax,subalgo_outer,subalgo_inner,subalgo_minmax.par_minmax);

        dmin_sp(dmin_sp < 0) = 0;
        dmin_sp(dmin_sp > 1) = 1;

        % Remove solutions returned more than once (if any)
        [~,idx] = unique(round(tol*dmin_sp),'rows','stable');
        dmin_sp=dmin_sp(idx,:);
        f_outer_sp = f_outer_sp(idx,:);
        for obj = 1:n_obj
            u_outer_sp{obj} = output_aux.u{obj}(idx,:);
        end
    end

    % outer loop with the archive chross.checks
    dmin_au = []; f_outer_au = [];
    if(algo_outer.par_au.use && (~alternate_outer || meth_outer(i_meth_outer) == 3))
        problem_min_d_au.par_objfun.u_record = u_record;

        if isfield (par_minmax, 'prob_reference_in_population') && rand()<par_minmax.prob_reference_in_population
            if n_obj > 1
                warning('reference in populaiton not implemented for MO')
            else
                algo_outer.par_au.initial_population = ref(:,1:n_d);
            end
        else
            algo_outer.par_au.initial_population = [];
            
        end

        [dmin_au, f_outer_au, ~ , output_aux_au] = algo_outer.optimise_au(problem_min_d_au,algo_outer.par_au);
        dmin_au(dmin_au < 0) = 0;
        dmin_au(dmin_au > 1) = 1;
        % Remove solutions returned more than once (if any)
        [~,idx] = unique(round(tol*dmin_au),'rows','stable');
        dmin_au = dmin_au(idx,:);
        f_outer_au = f_outer_au(idx,:);
     end

    % "OUTER" LOOP with the surrogate
    dmin_sa = []; f_outer_sa = [];
    if(algo_outer.par_sa.use&&(~alternate_outer || meth_outer(i_meth_outer) == 1))

        if isfield (par_minmax, 'prob_reference_in_population') && rand()<par_minmax.prob_reference_in_population
            if n_obj > 1
                warning('reference in populaiton not implemented for MO')
            else
                algo_outer.par_sa.initial_population = ref(:,1:n_d);
            end
        else
            algo_outer.par_sa.initial_population = [];
            
        end

        [dmin_sa, f_outer_sa, ~ , output_aux] = algo_outer.optimise_sa(problem_min_d_sa,algo_outer.par_sa);

        dmin_sa(dmin_sa < 0) = 0;
        dmin_sa(dmin_sa > 1) = 1;

        % Remove solutions returned more than once (if any)
        [~,idx] = unique(round(tol*dmin_sa),'rows','stable');
        dmin_sa=dmin_sa(idx,:);
        f_outer_sa = f_outer_sa(idx,:);
    end

    dmin = [dmin_sp; dmin_au; dmin_sa];
    f_outer = [f_outer_sp; f_outer_au; f_outer_sa];
    % Remove solutions returned more than once (if any)
    [~,idx] = unique(round(tol*dmin),'rows','stable');
    dmin=dmin(idx,:);
    f_outer = f_outer(idx,:);

    if (verbosity)
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
    end
    
    % "INNER" LOOP: Compute ue(i+1) = arg max f(dmin(i),u)
    record_aux = [];
    d_info_aux = [];
    f_inner = [];
    for i = 1:size(dmin,1)
        % even if you already know that d, you redo the maximisation bc you converged there
        di = dmin(i,:);
        problem_max_u.par_objfun.d = di;
        ufmax = nan(1,(n_u+1)*n_obj);

        for obj = 1:n_obj
            problem_max_u.par_objfun.objective = obj;

            % Give a member of the population by psi-surr.
            % Unless the guy comes from the S(d,u) in which case you priorize the u of S(d,u)
            algo_inner{obj}.par.initial_population= [];
            is_sp = false;
            loc_sp = 0;
            if ~isempty(dmin_sp)     
                [is_sp, loc_sp] = ismember(di,dmin_sp,'rows');
            end
            if is_sp
                u0 = u_outer_sp{obj}(loc_sp,:);
                algo_inner{obj}.par.initial_population = u0;

            elseif psi_surr(obj)
                [uf0,~] = problem_min_d_sa.par_objfun.surrogate.predictor(di, problem_min_d_sa.par_objfun.surrogate.model);
                u0 = uf0((obj-1)*n_u+1 : obj*n_u);
                u0(u0 < 0) = 0;
                u0(u0 > 1) = 1;
                algo_inner{obj}.par.initial_population= u0;
            end
                
            % Maximize subproblem to find umax = arg max f(dmin(i),u)
            [ umax, fmax , ~ , output_aux] = algo_inner{obj}.optimise(problem_max_u,algo_inner{obj}.par);
            % nfeval = nfeval + output_aux.nfeval;
            umax(umax < 0) = 0;
            umax(umax > 1) = 1;
            fmax = -sign_inner*fmax;
            
            % remember u and fmax
            ufmax((obj-1)*n_u+1 : obj*n_u) = umax;
            ufmax(n_u*n_obj+obj) = fmax;
            f_inner(i,obj) = fmax;
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
    
    % if ~isempty(record_aux)
    %     f_inner = record_aux(:,end-n_obj+1:end);
    % end
    if verbosity
        if (n_obj == 2)
            % figure(2)
            colors = 'ckbrmg';
            hold on
            plot(record_aux(:,end-1),record_aux(:,end),strcat(colors(mod(iter,length(colors))+1),'o'))
            drawnow
            % figure(1)    
        end
        f_inner
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
    if verbosity
        current_ref = ref;
        current_ref(:,1:n_d) = ref(:,1:n_d).*repmat(problem_minmax.ub_d'-problem_minmax.lb_d',[size(ref,1),1]) + repmat(problem_minmax.lb_d',[size(ref,1),1]);
        for obj = 1:n_obj
            map_info = problem_max_u.par_objfun.map_u_info{obj};
            for i = 1:size(ref,1)
                current_ref(i,n_d+(obj-1)*n_u+1:n_d+obj*n_u) = map_affine(ref(i,n_d+(obj-1)*n_u+1 : n_d+obj*n_u),map_info); %this can be easily vectorized
            end
        end
        current_ref
        nfevalglobal
        if (n_obj == 2)
            colors = 'ckbrmg';
            hold on
            plot(current_ref(:,end-1),current_ref(:,end),strcat(colors(mod(iter,length(colors))+1),'^'))
            title(num2str(nfevalglobal))
            drawnow
        end
        % figure(1)
        % plot(record(:,1),record(:,3),'.')
        % drawnow
        % figure(2)
        % plot(record(:,1),record(:,2),'.')
        % drawnow
    end
    
    % for alternation
    if par_minmax.alternate_outer
        if all(ismember(ref,last_ref,'rows'))
            i_meth_outer = i_meth_outer+1;
            iter_alternation = 0;
        else
            %i_meth_outer = 1;
            iter_alternation = iter_alternation+1;
            if(isfield(par_minmax,'max_iter_alternation') && iter_alternation >= par_minmax.max_iter_alternation)
                i_meth_outer = i_meth_outer+1;
                iter_alternation = 0;
            end
        end

        if i_meth_outer > length(meth_outer)
            i_meth_outer = 1;
        end

    end
    last_ref = ref;
    
    %format long
    convergence(end+1).iter = iter;
    convergence(end).nfeval = nfevalglobal;
    convergence(end).d = ref(:,1:n_d);
    convergence(end).u = ref(:,n_d+1:end-n_obj);
    convergence(end).f = ref(:,end-n_obj+1:end);
    stop = nfevalglobal >= nfevalmax || iter >= max_iter;
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
