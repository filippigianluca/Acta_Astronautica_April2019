function [d,fval,exitflag,output] = optimise_scal(problem_minmax, algo_outer, algo_inner, par_minmax)

verbose = true;

global nfevalglobal
optimise = @optimise_ausa_for_scal;
n_d = problem_minmax.dim_d;
n_u = problem_minmax.dim_u;
n_obj = problem_minmax.n_obj;
n_subproblems = par_minmax.n_subproblems;
sign_inner = problem_minmax.sign_inner; % 1 for minmax, -1 for minmin

d = nan(n_subproblems,n_d);
umax = nan(n_subproblems,n_u);
u_record=cell(1,n_obj);
for obj = 1:n_obj
    u_record{obj} = nan(n_subproblems,n_u);
end
fscal = nan(n_subproblems,1);
fval = nan(n_subproblems,n_obj);
nfeval = 0;
nfevalmax = floor((par_minmax.maxnfeval-n_subproblems*n_obj*algo_inner.par.nFeValMax)/par_minmax.n_subproblems); % AD HOC MPIDEA!

% solve single objective problems to get some ref points
frontspans = ones(1,n_obj);
utopia = zeros(1,n_obj);
lambdas = lambda_gen(n_subproblems,n_obj);

problem_max_u = build_metaproblem_macsminmax_inner(problem_minmax);

for subp = 1:n_obj
    problem_inner = build_metaproblem_scal_inner(problem_minmax,lambdas(subp,:), utopia, frontspans);
    [d(subp,:),fscal(subp,1),~,output_aux] = optimise(problem_minmax, algo_outer, algo_inner, par_minmax, problem_inner, nfevalmax)
    umax(subp,:) = output_aux.u;
    nfeval = nfeval + output_aux.nfeval;

    problem_max_u.par_objfun.d = d(subp,:);
    for obj = 1:n_obj
        problem_max_u.par_objfun.objective = obj;
        % Maximize subproblem to find umax = arg max f(dmin(subp),u)
        [ umax2, f_inner , ~ , output_aux] = algo_inner.optimise(problem_max_u,algo_inner.par);
        nfeval = nfeval + output_aux.nfeval;
        umax2(umax2 < 0) = 0;
        umax2(umax2 > 1) = 1;
        f_inner = -problem_minmax.sign_inner*f_inner;
        if subp == obj && sign_inner*f_inner < sign_inner*fscal(subp,1)
            fval(subp,obj) = fscal(subp,1);
            u_record{obj}(subp,:) = umax(subp,:);
        else
            fval(subp,obj) = f_inner;
            u_record{obj}(subp,:) = umax2;
        end
    end
end

% change utopia point
utopia_0 = fscal(1:n_obj,1)';
utopia_1 = diag(fval)';
frontspans = max(fval(1:n_obj,:),[],1)-min(fval(1:n_obj,:),[],1);
utopia = utopia_1-0.1*frontspans;
if verbose
    utopia_0
    utopia_1
    utopia
    if n_obj == 2
        hold on
        plot(utopia_0(1), utopia_0(2), 'c.');
        plot(utopia_1(1), utopia_1(2), 'c*');
        plot(utopia(1), utopia(2), 'c*');
        plot([utopia_0(1) utopia(1)], [utopia_0(2) utopia(2)], 'c--');
        plot(fval(1:n_obj,1),fval(1:n_obj,2),'r.')
        drawnow;
    end
end

% solve the rest of subproblems
if verbose
    frontspans
end

for subp = n_obj+1 : n_subproblems
    problem_inner = build_metaproblem_scal_inner(problem_minmax,lambdas(subp,:), utopia, frontspans);
    [d(subp,:),fscal(subp,1),~,output_aux] = optimise(problem_minmax, algo_outer, algo_inner, par_minmax, problem_inner, nfevalmax)
    umax(subp,:) = output_aux.u;
    nfeval = nfeval + output_aux.nfeval;

    problem_max_u.par_objfun.d = d(subp,:);
    for obj = 1:n_obj
        problem_max_u.par_objfun.objective = obj;
        % Maximize subproblem to find umax = arg max f(dmin(subp),u)
        [ umax2, f_inner , ~ , output_aux] = algo_inner.optimise(problem_max_u,algo_inner.par);
        nfeval = nfeval + output_aux.nfeval;
        umax2(umax2 < 0) = 0;
        umax2(umax2 > 1) = 1;
        f_inner = -problem_minmax.sign_inner*f_inner;
        fval(subp,obj) = f_inner;
        u_record{obj}(subp,:) = umax2;
        % if subp == obj && sign_inner*f_inner < sign_inner*f_scal(subp,1)
        %     fval(subp,obj) = f_scal(subp,1);
        % else
        %     fval(subp,obj) = f_inner;
        % end
    end
    if verbose
        lambdas(subp,:)
        fval
        if n_obj == 2
            plot(fval(subp,1),fval(subp,2),'r.')
            drawnow;
        end
    end
end

% if verbose
%     fval
%     if n_obj == 2
%         plot(fval(:,1),fval(:,2),'r.')
%         drawnow;
%     end
% end

d = d.*repmat(problem_minmax.ub_d'-problem_minmax.lb_d',[n_subproblems,1]) + repmat(problem_minmax.lb_d',[n_subproblems,1]);
u = cell(1,n_obj);
for obj = 1:n_obj
    map_info = problem_max_u.par_objfun.map_u_info{obj};
    for subp = 1:n_subproblems
        u{obj}(subp,:) = map_affine(u_record{obj}(subp,:),map_info); %this can be easily vectorized
    end
end

output.u = u;
output.nfeval = nfeval;
exitflag = 0;

return

function [d,fval,exitflag,output] = optimise_so_for_scal(problem_minmax, algo_outer, algo_inner, par_minmax, problem_inner,nfevalmax)

% Rename inputs
n_d = problem_minmax.dim_d;
n_u = problem_minmax.dim_u;
n_obj = 1;
sign_inner = problem_minmax.sign_inner; % 1 for minmax, -1 for minmin

% local search flags, ADD SANITY CHECKS AND STUFF HERE
lsflag_validation = par_minmax.local_search_flags.validation; % run either local search or func eval in validation
lsflag_inner = par_minmax.local_search_flags.inner;           % run either local search or func eval in SO problem
lsflag_outer = par_minmax.local_search_flags.outer;           % run either local search or func eval in MO problem


% Build metaproblems
problem_max_u = problem_inner;
problem_min_d = build_metaproblem_macsminmax_outer(problem_max_u,lsflag_outer);

% Random initial guess for d
n_d0 = par_minmax.n_d0;
d_0 = lhsu(zeros(1,n_d),ones(1,n_d),n_d0);

% Initialise archive d
d_record = [];
f_record = [];
% Now the fun begins
nfeval = 0;

% Initialise archive u
u_record = cell(1,n_obj);
for i = 1:size(d_0,1)
    problem_max_u.par_objfun.d = d_0(i,:);                                                 % tell metaproblem what d to fix
    for obj = 1:n_obj
        problem_max_u.par_objfun.objective = obj;                                     % tell metaproblem what objective to optimise
        [ umax, ~ , ~ , output_aux ] = algo_inner.optimise(problem_max_u,algo_inner.par); % optimise
        u_record{obj} = [u_record{obj}; umax];                                        % append result
        nfeval = nfeval + output_aux.nfeval;
    end
end

% Main loop
size_u_record = n_obj;
stop = false;
while ~stop
    problem_min_d.par_objfun.u_record = u_record;
    % "OUTER" LOOP: Compute dmin(i) = arg min {max f1(d,ue1), max f2(d,ue2)}
    [dmin, f_outer, ~ , output_aux] = algo_outer.optimise(problem_min_d,algo_outer.par);

    nfeval = nfeval + (output_aux.nfeval+(n_obj>1))*(~lsflag_outer + lsflag_outer*20*problem_minmax.dim_u)*size_u_record; 
    % NOTE: the +1 is empirical
    % NOTE: 50*dim_u is hardcoded nfevalmax in local search but is too conservative and overestimates nfeval hence 20.
    
    dmin(dmin < 0) = 0;
    dmin(dmin > 1) = 1;

    
    % Remove solutions archived more than once (if any)
    [~,idx] = unique(round(1e8*dmin),'rows');
    dmin=dmin(idx,:);
    f_outer = f_outer(idx,:);
    f_record_aux = f_outer;

    % "INNER" LOOP: Compute ue(i+1) = arg max f(dmin(i),u)
    for i = 1:size(dmin,1)
        problem_max_u.par_objfun.d = dmin(i,:);
        for obj = 1:n_obj
            problem_max_u.par_objfun.objective = obj;
            % Maximize subproblem to find umax = arg max f(dmin(i),u)
            [ umax, f_inner , ~ , output_aux] = algo_inner.optimise(problem_max_u,algo_inner.par);
            nfeval = nfeval + output_aux.nfeval;
            umax(umax < 0) = 0;
            umax(umax > 1) = 1;
            f_inner = -sign_inner*f_inner;
            
            % Archive best between umax and umax2 [in case IDEA failed]
            if sign_inner*f_inner > sign_inner*f_outer(i,obj)
                if lsflag_inner
                    u_cell_aux = cell(1,n_obj);
                    u_cell_aux{obj} = umax;
                    [fmax2,umax2,nfeval2] = u_validation(problem_max_u, dmin(i,:), u_cell_aux, lsflag_inner, obj);
                    nfeval = nfeval + nfeval2;

                    u_record{obj} = [u_record{obj}; umax2{obj}];
                    f_record_aux(i,obj) = fmax2(1,obj);
                else
                    u_record{obj} = [u_record{obj}; umax];
                    f_record_aux(i,obj) = f_inner;
                end 

                % If lsflag_inner = 0, the "else" step wastes evaluations and is
                % useless: It finds and archives a u that is already in the
                % archive, and is removed by the "unique" at line 96.
            else
                if lsflag_inner
                    % Evaluate subproblem to find umax2 = arg max f(dmin(i),ua)
                    [fmax2,umax2,nfeval2] = u_validation(problem_max_u, dmin(i,:), u_record, lsflag_inner, obj);
                    nfeval = nfeval + nfeval2;
                    u_record{obj} = [u_record{obj}; umax2{obj}];
                    f_record_aux(i,obj) = fmax2(1,obj);
                end
            end
        end
    end
    
    % Archive dmin
    d_record = [d_record; dmin];
    f_record = [f_record; f_record_aux];
    
    % Remove solutions archived more than once (if any)
    [~, idx] = unique(round(1e8*d_record),'rows');
    d_record = d_record(idx,:);
    f_record = f_record(idx,:); %not worrying too much by loss of validated solutions
    for obj = 1:n_obj
        [~, idx] = unique(round(1e8*u_record{obj}),'rows');
        u_record{obj} = u_record{obj}(idx,:);
    end
    
    size_u_record = sum(cellfun('size',u_record,1));
    nfeval_loop = nfevalmax - size_u_record*size(d_record,1)*(~lsflag_validation + lsflag_validation*20*problem_minmax.dim_u);
    % NOTE: 50*dim_u is hardcoded nfevalmax in local search but is too conservative and overestimates nfeval_val hence 20.

    
    stop = nfeval >= nfeval_loop;
end

%% archive cross check
% % Refine solutions that maximize the objectives. This step also associates a u vector to each d vector
% [ f_val, u_val_record, nfeval_val] = u_validation(problem_max_u, d_record, u_record, lsflag_validation, 1:n_obj);
% nfeval = nfeval + nfeval_val;

% %% Select non-dominated solutions and define outputs
% sel = dominance(fval,0) == 0;

% fval = f_val(sel,:);

%% more clever archive cross check
% even more clever would be to keep track of which u's have the d's been validated against already
stop = false;
checked = false(1,size(d_record,1));
sel = false(1,size(d_record,1));
u_val_record = cell(1,n_obj);
for obj = 1:n_obj
    u_val_record{obj}=nan(size(d_record,1),problem_minmax.dim_u);
end
nfeval_val = 0;
while ~stop
    sel = dominance(f_record,0) == 0;       % find non-dominated
    if(n_obj>1) 
        sel = sel';
    end
    
    tocheck = sel & ~checked;               % select only those that have not been checked
    stop = all(~tocheck);                   % stop if you have nothing to check
    % check those that have been selected
    d_tocheck = d_record(tocheck,:);
    [f_val_aux, u_val_record_aux, nfeval_aux] = u_validation(problem_max_u, d_tocheck, u_record, lsflag_validation, 1:n_obj);
    nfeval_val = nfeval_val + nfeval_aux;

    % update f_record and u_val_record
    f_record(tocheck,:) = f_val_aux;
    for obj = 1:n_obj
        u_val_record{obj}(tocheck,:) = u_val_record_aux{obj};
    end
    checked = checked | tocheck;           % update who has been checked
end

fval = f_record(sel,:);
nfeval = nfeval + nfeval_val;

%%
frontsize = size(fval,1);
d = d_record(sel,:);
for obj = 1:n_obj
    u_val_record{obj} = u_val_record{obj}(sel,:);
end

output.u = u_val_record{1};
output.nfeval = nfeval;
exitflag = 0;


% d = d_record(sel,:).*repmat(problem_minmax.ub_d'-problem_minmax.lb_d',[frontsize,1]) + repmat(problem_minmax.lb_d',[frontsize,1]);
% u = cell(1,n_obj);
% for obj = 1:n_obj
%     u_val_record{obj} = u_val_record{obj}(sel,:);
%     map_info = problem_max_u.par_objfun.map_u_info{obj};
%     for i = 1:frontsize
%         u{obj}(i,:) = map_affine(u_val_record{obj}(i,:),map_info); %this can be easily vectorized
%     end
% end

% output.u = u;
% output.nfeval = nfeval;
% exitflag = 0;

return

function [lambdas] = lambda_gen(n_subproblems,n_obj);

if n_obj == 2
    angles = linspace(0,pi/2,n_subproblems)';
    tans = tan(angles(2:end-1));
    lambdas = [eye(n_obj); 1./(1+tans) tans./(1+tans)];
else
    error('lambda_gen for n_obj>2 not yet implemented')
end


return

function [d,fval,exitflag,output] = optimise_ausa_for_scal(problem_minmax, algo_outer, algo_inner, par_minmax, problem_inner,nfevalmax)
verbose = false;
% Rename inputs
n_d = problem_minmax.dim_d;
n_u = problem_minmax.dim_u;
% n_obj = problem_minmax.n_obj;
n_obj = 1;
sign_inner = problem_minmax.sign_inner; % 1 for minmax, -1 for minmin

% nfevalmax = par_minmax.maxnfeval;

% local search flags, ADD SANITY CHECKS AND STUFF HERE
lsflag_validation = par_minmax.local_search_flags.validation; % run either local search or func eval in validation
lsflag_inner = par_minmax.local_search_flags.inner;           % run either local search or func eval in SO problem
lsflag_outer = par_minmax.local_search_flags.outer;


% Build metaproblems
% problem_max_u = build_metaproblem_macsminmax_inner(problem_minmax);
problem_max_u = problem_inner;
problem_min_d_au = build_metaproblem_macsminmax_outer(problem_max_u,lsflag_outer);
problem_min_d_sa = build_metaproblem_mo_sa_outer(problem_minmax);

% Random initial guess for d
n_d0 = par_minmax.n_d0;
d_0 = lhsu(zeros(1,n_d),ones(1,n_d),n_d0);

% Initialise archive d
d_record = d_0;
f_record = nan(n_d0,n_obj);
% Now the fun begins
nfeval = 0;

% Initialise archive u
u_record = cell(1,n_obj);
for i = 1:size(d_0,1)
    problem_max_u.par_objfun.d = d_0(i,:);                                                 % tell metaproblem what d to fix
    for obj = 1:n_obj
        problem_max_u.par_objfun.objective = obj;                                     % tell metaproblem what objective to optimise
        [ umax, f_aux, ~ , output_aux ] = algo_inner.optimise(problem_max_u,algo_inner.par); % optimise
        u_record{obj} = [u_record{obj}; umax];                                        % append result
        f_record(i,obj) = -sign_inner*f_aux;
        nfeval = nfeval + output_aux.nfeval;
    end
end

% Remove solutions archived more than once (if any)
[~, idx] = unique(round(1e8*d_record),'rows');
d_record = d_record(idx,:);
f_record = f_record(idx,:); %not worrying too much by loss of validated solutions
for obj = 1:n_obj
    [~, idx] = unique(round(1e8*u_record{obj}),'rows');
    u_record{obj} = u_record{obj}(idx,:);
end

% Main loop
iter = 0;
size_u_record = n_obj;
stop = false;
while ~stop
    iter=iter+1;

    % "OUTER" LOOP: Compute dmin(i) = arg min {max f1(d,ue1), max f2(d,ue2)}
   
    % Au minimisation
    problem_min_d_au.par_objfun.u_record = u_record;
    dmin_au = []; f_outer_au = [];
    if(algo_outer.par_au.use)
        [dmin_au, f_outer_au, ~ , output_aux_au] = algo_outer.optimise_au(problem_min_d_au,algo_outer.par_au);
        nfeval = nfeval + (output_aux_au.nfeval+(n_obj>1))*(~lsflag_outer + lsflag_outer*20*problem_minmax.dim_u)*size_u_record; 
        % NOTE: the +1 is empirical
        % NOTE: 50*dim_u is hardcoded nfevalmax in local search but is too conservative and overestimates nfeval hence 20.
        dmin_au(dmin_au < 0) = 0;
        dmin_au(dmin_au > 1) = 1;
        [~,idx] = unique(round(1e8*dmin_au),'rows');
        dmin_au=dmin_au(idx,:);
        f_outer_au = f_outer_au(idx,:);
    end
    
    dmin_sa = []; f_outer_sa = [];
    if(algo_outer.par_sa.use)
        % surrogate-based minimisation. The surrogate is not trained after the Au minimisation to avoid non-maximised solutions in the surrogate.
        problem_min_d_sa.par_objfun.surrogate.model = problem_min_d_sa.par_objfun.surrogate.training(d_record,f_record,problem_min_d_sa.par_objfun.surrogate);
        % problem_min_d.par_objfun.surrogate.model = problem_min_d.par_objfun.surrogate.training(dmin,f_record_aux,problem_min_d.par_objfun.surrogate);
        [dmin_sa, f_outer_sa, ~ , output_aux] = algo_outer.optimise_sa(problem_min_d_sa,algo_outer.par_sa);
        dmin_sa(dmin_sa < 0) = 0;
        dmin_sa(dmin_sa > 1) = 1;
        % Remove solutions archived more than once (if any)
        [~,idx] = unique(round(1e8*dmin_sa),'rows');
        dmin_sa=dmin_sa(idx,:);
        f_outer_sa = f_outer_sa(idx,:);

    end
    
    f_record_aux = [f_outer_sa; f_outer_au];
    
    if(verbose)
        if (n_obj == 2)
            % figure(2)
            colors = 'ckbrmg';
            hold on
            if(algo_outer.par_sa.use)
                plot(f_outer_sa(:,1),f_outer_sa(:,2),strcat(colors(mod(iter,length(colors))+1),'.'))
            end
            if(algo_outer.par_au.use)
                plot(f_outer_au(:,1),f_outer_au(:,2),strcat(colors(mod(iter,length(colors))+1),'s'))
            end
            drawnow
            % figure(1)    
        end
        f_record_aux
    end


    % "INNER" LOOP: Compute ue(i+1) = arg max f(dmin(i),u)
    for i = 1:size(dmin_sa,1)
        problem_max_u.par_objfun.d = dmin_sa(i,:);
        for obj = 1:n_obj
            problem_max_u.par_objfun.objective = obj;
            % Maximize subproblem to find umax = arg max f(dmin(i),u)
            [ umax, f_inner , ~ , output_aux] = algo_inner.optimise(problem_max_u,algo_inner.par);
            nfeval = nfeval + output_aux.nfeval;
            umax(umax < 0) = 0;
            umax(umax > 1) = 1;
            f_inner = -sign_inner*f_inner;
            
            % ALWAYS ARCHIVE
            u_record{obj} = [u_record{obj}; umax];
            f_record_aux(i,obj) = f_inner;
            
%             % Archive best between umax and umax2 [in case IDEA failed]
%             if sign_inner*f_inner > sign_inner*f_outer(i,obj)
%                 if lsflag_inner
%                     u_cell_aux = cell(1,n_obj);
%                     u_cell_aux{obj} = umax;
%                     [fmax2,umax2,nfeval2] = u_validation(problem_max_u, dmin(i,:), u_cell_aux, lsflag_inner, obj);
%                     nfeval = nfeval + nfeval2;
% 
%                     u_record{obj} = [u_record{obj}; umax2{obj}];
%                     f_record_aux(i,obj) = fmax2(1,obj);
%                 else
%                     u_record{obj} = [u_record{obj}; umax];
%                     f_record_aux(i,obj) = f_inner;
%                 end 
% 
%                 % If lsflag_inner = 0, the "else" step wastes evaluations and is
%                 % useless: It finds and archives a u that is already in the
%                 % archive, and is removed by the "unique" at line 96.
%             else
%                 if lsflag_inner
%                     % Evaluate subproblem to find umax2 = arg max f(dmin(i),ua)
%                     [fmax2,umax2,nfeval2] = u_validation(problem_max_u, dmin(i,:), u_record, lsflag_inner, obj);
%                     nfeval = nfeval + nfeval2;
%                     u_record{obj} = [u_record{obj}; umax2{obj}];
%                     f_record_aux(i,obj) = fmax2(1,obj);
%                 end
%             end
        end
    end
    
    % "INNER" LOOP AU: Compute ue(i+1) = arg max f(dmin(i),u)
    % u_record_iter = cell(1,n_obj);
    for i = 1:size(dmin_au,1)
        problem_max_u.par_objfun.d = dmin_au(i,:);
        for obj = 1:n_obj
            problem_max_u.par_objfun.objective = obj;
            % Maximize subproblem to find umax = arg max f(dmin(i),u)
            [ umax, f_inner , ~ , output_aux] = algo_inner.optimise(problem_max_u,algo_inner.par);
            nfeval = nfeval + output_aux.nfeval;
            umax(umax < 0) = 0;
            umax(umax > 1) = 1;
            f_inner = -sign_inner*f_inner;
            
            % Archive best between umax and umax2 [in case IDEA failed]
            if sign_inner*f_inner > sign_inner*f_outer_au(i,obj)
                if lsflag_inner(obj)
                    u_cell_aux = cell(1,n_obj);
                    u_cell_aux{obj} = umax;
                    [fmax2,umax2,nfeval2] = u_validation(problem_max_u, dmin_au(i,:), u_cell_aux, lsflag_inner(obj), obj);
                    nfeval = nfeval + nfeval2;

                    u_record{obj} = [u_record{obj}; umax2{obj}];
                    f_record_aux(size(dmin_sa,1)+i,obj) = fmax2(1,obj);
                else
                    % u_record_iter{obj} = [u_record_iter{obj}; umax];
                    u_record{obj} = [u_record{obj}; umax];
                    f_record_aux(size(dmin_sa,1)+i,obj) = f_inner;
                end 

                % If lsflag_inner = 0, the "else" step wastes evaluations and is
                % useless: It finds and archives a u that is already in the
                % archive, and is removed by the "unique" later
            else
                if lsflag_inner(obj)
                    % Evaluate subproblem to find umax2 = arg max f(dmin(i),ua)
                    [fmax2,umax2,nfeval2] = u_validation(problem_max_u, dmin_au(i,:), u_record, lsflag_inner(obj), obj);
                    nfeval = nfeval + nfeval2;
                    u_record{obj} = [u_record{obj}; umax2{obj}];
                    f_record_aux(size(dmin_sa,1)+i,obj) = fmax2(1,obj);
                % else
                %     [~,umax2,~] = u_validation(problem_max_u, dmin_au(i,:), u_record, false, obj);
                %     u_record_iter{obj} = [u_record_iter{obj}; umax2{obj}];
                end
            end
        end
    end
    
    if(verbose)
        if (n_obj == 2)
            % figure(2)
            colors = 'ckbrmg';
            hold on
            plot(f_record_aux(:,1),f_record_aux(:,2),strcat(colors(mod(iter,length(colors))+1),'o'))
            drawnow
            % figure(1)    
        end
        f_record_aux
    end

    % Archive dmin
    d_record = [d_record; dmin_sa ; dmin_au];
    f_record = [f_record; f_record_aux];
    
    % Remove solutions archived more than once (if any)
    [~, idx] = unique(round(1e8*d_record),'rows');
    d_record = d_record(idx,:);
    f_record = f_record(idx,:); %not worrying too much by loss of validated solutions
    for obj = 1:n_obj
        [~, idx] = unique(round(1e8*u_record{obj}),'rows');
        u_record{obj} = u_record{obj}(idx,:);
    end

    [f_record, u_val_record_aux, nfeval_aux] = u_validation(problem_max_u, d_record, u_record, false, 1:n_obj);
    nfeval = nfeval + nfeval_aux;

    % % remove from u_record the solutions that do not maximise any d in d_record...?
    % u_record = u_val_record_aux;
    % for obj = 1:n_obj
    %     [~, idx] = unique(round(1e8*u_record{obj}),'rows');
    %     u_record{obj} = u_record{obj}(idx,:);
    % end

    size_u_record = sum(cellfun('size',u_record,1));
    nfeval_loop = nfevalmax - size_u_record*size(d_record,1)*(lsflag_validation*20*problem_minmax.dim_u);
    % NOTE: 50*dim_u is hardcoded nfevalmax in local search but is too conservative and overestimates nfeval_val hence 20.

    stop = nfeval >= nfeval_loop;
end

%% archive cross check
% % Refine solutions that maximize the objectives. This step also associates a u vector to each d vector
% [ f_val, u_val_record, nfeval_val] = u_validation(problem_max_u, d_record, u_record, lsflag_validation, 1:n_obj);
% nfeval = nfeval + nfeval_val;

% %% Select non-dominated solutions and define outputs
% sel = dominance(fval,0) == 0;

% fval = f_val(sel,:);

%% more clever archive cross check
% even more clever would be to keep track of which u's have the d's been validated against already
sel = false(1,size(d_record,1));
u_val_record = cell(1,n_obj);
fval = [];
if (lsflag_validation)
    stop = false;
    checked = false(1,size(d_record,1));
    for obj = 1:n_obj
        u_val_record{obj}=nan(size(d_record,1),problem_minmax.dim_u);
    end
    nfeval_val = 0;
    while ~stop
        sel = dominance(f_record,0) == 0;       % find non-dominated
        if(n_obj>1) 
            sel = sel';
        end
        
        tocheck = sel & ~checked;               % select only those that have not been checked
        stop = all(~tocheck);                   % stop if you have nothing to check
        % check those that have been selected
        d_tocheck = d_record(tocheck,:);
        [f_val_aux, u_val_record_aux, nfeval_aux] = u_validation(problem_max_u, d_tocheck, u_record, lsflag_validation, 1:n_obj);
        nfeval_val = nfeval_val + nfeval_aux;

        % update f_record and u_val_record
        f_record(tocheck,:) = f_val_aux;
        for obj = 1:n_obj
            u_val_record{obj}(tocheck,:) = u_val_record_aux{obj};
        end
        checked = checked | tocheck;           % update who has been checked
    end

    fval = f_record(sel,:);
    nfeval = nfeval + nfeval_val;
else
    sel = dominance(f_record,0) == 0;       % find non-dominated
    fval = f_record(sel,:);
    u_val_record = u_val_record_aux;
end


%%
frontsize = size(fval,1);
d = d_record(sel,:);
for obj = 1:n_obj
    u_val_record{obj} = u_val_record{obj}(sel,:);
end

output.u = u_val_record{1};
output.nfeval = nfeval;
exitflag = 0;


% d = d_record(sel,:).*repmat(problem_minmax.ub_d'-problem_minmax.lb_d',[frontsize,1]) + repmat(problem_minmax.lb_d',[frontsize,1]);
% u = cell(1,n_obj);
% for obj = 1:n_obj
%     u_val_record{obj} = u_val_record{obj}(sel,:);
%     map_info = problem_max_u.par_objfun.map_u_info{obj};
%     for i = 1:frontsize
%         u{obj}(i,:) = map_affine(u_val_record{obj}(i,:),map_info); %this can be easily vectorized
%     end
% end

% output.u = u;
% output.nfeval = nfeval;
% exitflag = 0;




return