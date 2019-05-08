% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [d,fval,exitflag,output] = optimise_so_run_unconstr_nu400_20000inout(problem_minmax, algo_outer, algo_inner, par_minmax)

N_reduce_record_constraint = 0;

% Rename inputs
n_d = problem_minmax.dim_d;
n_u = problem_minmax.dim_u;
n_obj = problem_minmax.n_obj;


problem_minmax.sign_inner = 1;
sign_inner = problem_minmax.sign_inner; % 1 for minmax, -1 for minmin



nfevalmax = par_minmax.maxnfeval;

N_u_mpaidea = 2;
N_txt_mpaidea_inner  = 1;
N_txt_mpaidea_constr = 1;
N_txt_mpaidea_outer  = 1;


EPS_constraint = 0;  % realmin('single');

% doesn't accept the relaxation
flag_accept_relaxation = 0;




% open txt files
[str_archive, str_archive_inner, str_archive_outer] = open_txt_optimise_so(problem_minmax, algo_inner);



% local search flags, ADD SANITY CHECKS AND STUFF HERE
lsflag_validation = par_minmax.local_search_flags.validation; % run either local search or func eval in validation
lsflag_inner = par_minmax.local_search_flags.inner;           % run either local search or func eval in SO problem
lsflag_outer = par_minmax.local_search_flags.outer;           % run either local search or func eval in MO problem

algo_inner.par.sign_inner=sign_inner;
algo_outer.par.sign_inner=sign_inner;



% Build metaproblems
problem_max_u = build_metaproblem_macsminmax_inner(problem_minmax);
problem_min_d = build_metaproblem_macsminmax_outer(problem_max_u,lsflag_outer);


% problem_max_u.par_objfun.problem_par_objfun{n_obj}.fix=problem_minmax.fix;       % FIXED PARAMETERS

% Random initial guess for d
n_d0 = par_minmax.n_d0;
d_0 = lhsu(zeros(1,n_d),ones(1,n_d),n_d0);

% Initialise archive d
d_record = [];
f_record = [];

f_constraint_record = [];


nfeval = 0; % number of objective-function evaluation
nceval = 0; % number of constraint-function evaluation
n_eval = 0; % number of function evaluation: I have to decide whether (nceval + nfeval) or other ...
% at the moment:
%                n_eval = nfeval;


%----------------------------------------------------------------------
% INNER LOOP: Initialise archive u
%----------------------------------------------------------------------

u_record_maxC = cell(1, n_obj);
u_record      = cell(1, n_obj);
ARCHIVE{1}    = [{'d_min'}, {'u_max'}, {'u_max_C'}, {'f_inner'}, {'violation_C'}, {'max_violation_C'}, {'N_feval'}];
ARCHIVE_normalised{1}    = [{'d_min'}, {'u_max'}, {'u_max_C'}, {'f_inner'}, {'violation_C'}, {'max_violation_C'}, {'N_feval'}];
ARCHIVE_C{1}  = [{'d_min'}, {'u_max_C'}, {'f_max_C'}];
ARCHIVE_INNER = cell(1, n_obj);
ARCHIVE_OUTER = cell(1, n_obj);
umax_constraint = [];
violation_kk    = [];
violation_d_kk  = [];
violation_u_kk  = [];
f_inner_constraint = [];

for i = 1:size(d_0,1)
    problem_max_u.par_objfun.d = d_0(i,:);                                                % tell metaproblem what d to fix
    for obj = 1:n_obj
        
        
        algo_inner = save_MP_AIDEA_inner_openfile(problem_minmax, N_txt_mpaidea_inner, algo_inner);  % save population MPAIDEA
        
        problem_max_u.par_objfun.objective = obj;                                         % tell metaproblem what objective to optimise
        [ umax_all, fmax1_all , ~ , output_aux ] = algo_inner.optimise(problem_max_u,algo_inner.par); % optimise
        
        map_info = problem_max_u.par_objfun.map_u_info{obj};
        
        % for the constrained problem
        if ~isempty(problem_max_u.constraint) || problem_max_u.par_objfun.obj_constr
            for kk = 1:size(umax_all, 1)
                d = d_0(i,:).*(problem_minmax.ub_d'-problem_minmax.lb_d') + (problem_minmax.lb_d');
                u = map_affine(umax_all(kk,:), map_info);
                if ~isempty(problem_max_u.constraint)
                    violation(kk) = problem_minmax.constraint{obj}(d, u, problem_max_u.par_objfun.problem_par_objfun{obj});
                else
                    [out_obj_constr] = problem_minmax.objfun{obj}(d, u, problem_max_u.par_objfun.problem_par_objfun{obj});
                    fmax1_all(kk) = - out_obj_constr.f;
                    violation(kk) = out_obj_constr.c;
                end
                if violation(kk) > EPS_constraint
                    
                    if ~isempty(problem_max_u.constraint)
                        fmax1_all(kk) = -problem_minmax.objfun{obj}(d, u, problem_max_u.par_objfun.problem_par_objfun{obj});
                    end
                    
                    nfeval = nfeval + 1;
                end
                
            end
            
            if all(violation > EPS_constraint) % to change: if violations are = --> choose the min arg min(-f)
                [violation_kk, n_max] = min(violation);
                umax = umax_all(n_max, :);
                fmax1 = fmax1_all(n_max);
                % %---------------------------------------------------------------
                %                     [~, n_maxC] = max(violation);
                %                     umaxC = umax_all(n_maxC, :);
                %                     fmaxC = fmax1_all(n_maxC);
                %                     f_constraint_record = [f_constraint_record; fmaxC];
                %                     u_record_maxC{obj} = [u_record_maxC{obj}; umaxC];
                % %------------------------------------------------------------
                
            else
                pos_leq_eps = find(violation <= 0);
                fmax1_all_feasible = fmax1_all(pos_leq_eps);
                [fmax1, n_max] = min(fmax1_all_feasible);
                umax = umax_all((fmax1_all==fmax1), :);
                umax = umax(1,:);
                violation_kk = violation(pos_leq_eps(n_max));
            end
            
            ARCHIVE_INNER{obj}  = [ARCHIVE_INNER{obj}; -fmax1 violation_kk];
            
            
            % for the un-constrained problem
        else
            [fmax1, n_max] = min(fmax1_all);
            umax = umax_all(n_max, :);
            
            ARCHIVE_INNER{obj}  = [ARCHIVE_INNER{obj}; -fmax1 ];
        end
        
        
        
        % save population MPAIDEA
        [~] = save_MP_AIDEA_closefile(algo_inner);
        
        N_txt_mpaidea_inner = N_txt_mpaidea_inner + 1;
        
        
        
        u_record{obj} = [u_record{obj}; umax];                                            % append result
        [maxfmax1, ~] = max(-fmax1);
        nfeval = nfeval + output_aux.nfeval;
        
        
        
        
        
        
        
        
        %----------------------------------------------------------------------
        % CONSTRAINTS
        if ~isempty(problem_max_u.constraint) || problem_minmax.obj_constr
            
            % save population MPAIDEA
            algo_inner = save_MP_AIDEA_constraint_openfile(problem_minmax, N_txt_mpaidea_constr, algo_inner);
            
            [ umax_constraint, f_inner_constraint , ~ , output_aux] = optimise_constraint(problem_max_u,algo_inner.par);
            % "optimise_constraint" change the sign of C:
            % max violation of C = min(- C)
            
            nceval = nceval + output_aux.nfeval;
            
            
            umax_constraint(umax_constraint < 0) = 0;
            umax_constraint(umax_constraint > 1) = 1;
            f_inner_constraint = sign_inner*f_inner_constraint;
            
            
            % add to u record the u that maximise the constraint violation
            if any(f_inner_constraint > EPS_constraint)
                f_constraint_record = [f_constraint_record; f_inner_constraint];
                u_record_maxC{obj} = [u_record_maxC{obj}; umax_constraint];
            end
            
            % save population MPAIDEA
            [~] = save_MP_AIDEA_closefile(algo_inner);
            N_txt_mpaidea_constr = N_txt_mpaidea_constr + 1;
            
        end
        %----------------------------------------------------------------------
        % save exact valoues of design d and uncertaint u
        d_0_true  = d_0(i,:).*(problem_minmax.ub_d'-problem_minmax.lb_d') + (problem_minmax.lb_d');
        umax_true = map_affine(umax, map_info);
        
        
        [ARCHIVE, ARCHIVE_normalised] = archive_archive_normalised(ARCHIVE, ARCHIVE_normalised, problem_max_u, map_info, ...
            umax_constraint, d_0_true, umax_true, maxfmax1, violation_kk, ...
            f_inner_constraint, nfeval, d_0, umax, obj);
        
        
        if algo_inner.save_archive
            
            fprintf(fileID_archive, str_archive, ARCHIVE{obj}(end,:));
            fprintf(fileID_archive_inner, str_archive_inner, ARCHIVE_INNER{obj}(end,:));
        end
        
    end
end
%----------------------------------------------------------------------



global loop
loop = 0;

%--------------------------------------------------------------------------
% MAIN LOOP
%--------------------------------------------------------------------------


stop = false;
flag_relaxation = 0;

while ~stop
    clearvars u_record_tot
    clearvars violation_d
    clearvars violation_u
    clearvars fmax1_all_feasible
    clearvars max_violation_d
    
    
    loop = loop + 1;
    
    
    for obj = 1:n_obj
        [~, idx] = unique(round(1e8*u_record{obj}),'rows');
        u_record{obj} = u_record{obj}(idx,:);
        
        [~, idx] = unique(round(1e8*u_record_maxC{obj}),'rows');
        u_record_maxC{obj} = u_record_maxC{obj}(idx,:);
        
        f_constraint_record = f_constraint_record(idx);
    end
    
    u_record_tot = {[u_record{obj}; u_record_maxC{obj}]};
    
    for obj = 1:n_obj
        [~, idx] = unique(round(1e8*u_record_tot{obj}),'rows');
        u_record_tot{obj} = u_record_tot{obj}(idx,:);
    end
    
    size_u_record = sum(cellfun('size',u_record_tot,1));
    
    
    problem_min_d.par_objfun.u_record = u_record_tot;
    
    
    % FIXED PARAMETERS
%     problem_min_d.par_objfun.problem_fix_d.par_objfun.problem_par_objfun{n_obj}.fix = problem_minmax.fix;
%     problem_min_d.par_objfun.problem_par_objfun{n_obj}.fix  = problem_minmax.fix;
%     problem_min_d.par_objfun.problem_par_objfun{n_obj}.ub_u = problem_minmax.ub_u{1};
 
    
    
    %--------------------------------------------------------------------------
    % "OUTER" LOOP: Compute dmin(i) = arg min {max f1(d,ue1), max f2(d,ue2)}
    
    % save population MPAIDEA
    algo_outer = save_MP_AIDEA_outer_loop_openfile(problem_minmax, N_txt_mpaidea_outer, algo_outer);
    
    [dmin_all, f_outer_all, ~ , output_aux] = algo_outer.optimise(problem_min_d,algo_outer.par);
    
    
    if ~isempty(problem_max_u.constraint) || problem_max_u.par_objfun.obj_constr
        
        max_violation_d = [];
        for kk=1:size(dmin_all, 1)
            d = dmin_all(kk,:).*(problem_minmax.ub_d'-problem_minmax.lb_d') + (problem_minmax.lb_d');
            map_info = problem_max_u.par_objfun.map_u_info{obj};
            violation_d = zeros(1,size(u_record_tot{obj},1));
            fmax_all    = [];
            violation_d = [];
            for ukk=1:size(u_record_tot{obj},1)
                u = map_affine(u_record_tot{obj}(ukk,:), map_info);
                
                
                if ~isempty(problem_max_u.constraint)
                    violation_d(ukk) = problem_minmax.constraint{obj}(d, u, problem_max_u.par_objfun.problem_par_objfun{obj});
                else
                    [out_obj_constr] = problem_minmax.objfun{obj}(d, u, problem_max_u.par_objfun.problem_par_objfun{obj});
                    fmax_all(ukk) =   out_obj_constr.f;
                    violation_d(ukk) = out_obj_constr.c;
                end
                
                
                nceval = nceval + 1;
            end
            [max_viol, pos_max_viol] = max(violation_d);
            max_violation_d(kk) = max_viol;
            if max_violation_d(kk) > EPS_constraint
                
                    u = map_affine(u_record_tot{obj}(pos_max_viol,:), map_info);
                    if ~isempty(problem_max_u.constraint)
                        f_outer_all(kk) = problem_minmax.objfun{obj}(d, u, problem_max_u.par_objfun.problem_par_objfun{obj});
                    else
                        f_outer_all(kk) = max(fmax_all);
                    end
                
                
                nfeval = nfeval + 1;
            end
        end
        
        if all(max_violation_d > EPS_constraint)
            [violation_d_kk, n_min] = min(max_violation_d);
            dmin = dmin_all(n_min, :);
            f_outer = f_outer_all(n_min);
            % %---------------------------------------------------------------
            %                 [~, n_minC] = max(max_violation_d);
            %                 uminC = dmin_all(n_minC, :);
            %                 fminC = f_outer_all(n_minC);
            %                 f_constraint_record = [f_constraint_record; fminC];
            %                 u_record_maxC{obj} = [u_record_maxC{obj}; uminC];
            % %------------------------------------------------------------
        else
            pos_leq_eps = find(max_violation_d <= EPS_constraint);
            fmax1_all_feasible = f_outer_all(pos_leq_eps);
            [f_outer, n_min] = min(fmax1_all_feasible);
            dmin = dmin_all((f_outer_all==f_outer), :);
            dmin = dmin(1,:);
            violation_d_kk = max_violation_d(pos_leq_eps(n_min));
        end
        
        ARCHIVE_OUTER{obj}  = [ARCHIVE_OUTER{obj}; f_outer violation_d_kk];
        
        
    else    % if no constraint
        [f_outer, n_min] = min(f_outer_all);
        dmin = dmin_all(n_min, :);
        
        ARCHIVE_OUTER{obj}  = [ARCHIVE_OUTER{obj}; f_outer ];
    end
    
    
    
    
    
    % save population MPAIDEA
    [~] = save_MP_AIDEA_closefile(algo_outer);
    
    N_txt_mpaidea_outer = N_txt_mpaidea_outer + 1;
    
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
    
    %         global fout
    %         fout = [fout f_outer];
    
    
    
    %--------------------------------------------------------------------------
    % "INNER" LOOP: Compute ue(i+1) = arg max f(dmin(i),u)
    for i = 1:size(dmin,1)
        problem_max_u.par_objfun.d = dmin(i,:);
        for obj = 1:n_obj
            
            % 1) evaluate the feasible  max(F(u))
            %
            
            
            % save population MPAIDEA
            algo_inner = save_MP_AIDEA_inner_loop_openfile(problem_minmax, N_txt_mpaidea_inner, algo_inner);
            
            
            problem_max_u.par_objfun.objective = obj;
            % Maximize subproblem to find umax = arg max f(dmin(i),u)
            [ umax_all, f_inner_all , ~ , output_aux] = algo_inner.optimise(problem_max_u,algo_inner.par);
            
            
            if ~isempty(problem_max_u.constraint) || problem_max_u.par_objfun.obj_constr
                for kk=1:size(umax_all,1)
                    d = dmin.*(problem_minmax.ub_d'-problem_minmax.lb_d') + (problem_minmax.lb_d');
                    map_info = problem_max_u.par_objfun.map_u_info{obj};
                    u = map_affine(umax_all(kk,:), map_info);
                    
                    
                    
                if ~isempty(problem_max_u.constraint)
                    violation_u(kk) = problem_minmax.constraint{obj}(d, u, problem_max_u.par_objfun.problem_par_objfun{obj});
                else
                    [out_obj_constr] = problem_minmax.objfun{obj}(d, u, problem_max_u.par_objfun.problem_par_objfun{obj});
                    f_inner_all(kk) = - out_obj_constr.f;
                    violation_u(kk) = out_obj_constr.c;
                end              
                    nceval = nceval + 1;
                    
                    if ~isempty(problem_max_u.constraint)
                        if violation_u(kk) ~= 0
                            f_inner_all(kk) = -problem_minmax.objfun{obj}(d, u, problem_max_u.par_objfun.problem_par_objfun{obj});

                            nfeval = nfeval + 1;
                        end
                    end
                end
                
                
                if all(violation_u > EPS_constraint)
                    [violation_u_kk, n_max] = min(violation_u);
                    umax = umax_all(n_max, :);
                    f_inner = f_inner_all(n_max);
                    % %---------------------------------------------------------------
                    %                         [~, n_maxC] = max(violation_u);
                    %                         umaxC = umax_all(n_maxC, :);
                    %                         fmaxC = fmax1_all(n_maxC);
                    %                         f_constraint_record = [f_constraint_record; fmaxC];
                    %                         u_record_maxC{obj} = [u_record_maxC{obj}; umaxC];
                    % %------------------------------------------------------------
                else
                    pos_leq_eps = find(violation_u <= EPS_constraint);
                    fmax1_all_feasible = f_inner_all(pos_leq_eps);
                    [f_inner, n_max] = min(fmax1_all_feasible);
                    umax = umax_all((f_inner_all==f_inner), :);
                    umax = umax(1, :);
                    violation_u_kk = violation_u(pos_leq_eps(n_max));
                end
                
                
                %                 [f_inner, n_max] = min(f_inner_all);
                %                 umax = umax_all(n_max, :);
                
                
                ARCHIVE_INNER{obj}  = [ARCHIVE_INNER{obj}; -f_inner violation_u_kk];
                
                
            else
                [f_inner, n_max] = min(f_inner_all);
                umax = umax_all(n_max, :);
                
                ARCHIVE_INNER{obj}  = [ARCHIVE_INNER{obj}; -f_inner ];
            end
            
            
            
            %%save population MPAIDEA
            [~] = save_MP_AIDEA_closefile(algo_inner);
            
            N_txt_mpaidea_inner = N_txt_mpaidea_inner + 1;
            
            
            nfeval = nfeval + output_aux.nfeval;
            umax(umax < 0) = 0;
            umax(umax > 1) = 1;
            f_inner = -sign_inner*f_inner;
            
            
            
            
            
            %----------------------------------------------------------
            % CONSTRAINTS
            if ~isempty(problem_max_u.constraint) || problem_minmax.obj_constr
                % 2) evaluate the max violation of the constraint
                
                
                % save population MPAIDEA
                algo_inner = save_MP_AIDEA_constraint_loop_openfile(problem_minmax, N_txt_mpaidea_constr, algo_inner);
                
                
                [ umax_constraint, f_inner_constraint , ~ , output_aux] = optimise_constraint(problem_max_u,algo_inner.par);
                % "optimise_constraint" change the sign of C:
                % max violation of C = min(- C)
                
                %%save population MPAIDEA
                [~] = save_MP_AIDEA_closefile(algo_inner);
                N_txt_mpaidea_constr = N_txt_mpaidea_constr + 1;
                
                nfeval = nfeval + output_aux.nfeval;
                
                umax_constraint(umax_constraint < 0) = 0;
                umax_constraint(umax_constraint > 1) = 1;
                f_inner_constraint = sign_inner*f_inner_constraint;
                
                
                % chose the unfeasible u if the constraint is not respected in all the domain.
                
                % constraint relaxation
                N_eps_relaxation = 4;
                if ~isempty(f_constraint_record)
                    epsilon_min = min(min(f_constraint_record));
                    epsilon_max = max(max(f_constraint_record));
                else
                    epsilon_min = 0;
                    epsilon_max = 0;
                end
                
                OneOtherLoop = (algo_outer.par.nFeValMax*(size_u_record+1*(f_inner_constraint>0)) + 2*algo_inner.par.nFeValMax);
                
                Delta_f_inner_min = epsilon_min*(nfevalmax - nfeval < N_eps_relaxation*OneOtherLoop);
                Delta_f_inner_max = epsilon_max*(nfevalmax - nfeval < N_eps_relaxation*OneOtherLoop);
                %                         epsilon*(nfevalmax - nfeval < 3*OneOtherLoop) + ...
                %                         epsilon*(nfevalmax - nfeval < 2*OneOtherLoop) + ...
                %                         epsilon*(nfevalmax - nfeval < 1*OneOtherLoop);
                
                %                     Archive_matrix_constraint = cell2mat(ARCHIVE{obj}(2:end,:));
                
                if  ( ( nfevalmax - nfeval > N_eps_relaxation*OneOtherLoop || ...
                        any(f_constraint_record <= 0) ) && ...
                        flag_relaxation == 0 ) || flag_accept_relaxation == 0
                    
                    
                    % f_inner = f_inner_constraint;
                    % if the constraint is violated the corresponding u
                    % is recorded
                    
                    f_constraint_record = [f_constraint_record; f_inner_constraint];
                    u_record_maxC{obj} = [u_record_maxC{obj}; umax_constraint];
                    
                    
                elseif   ( ( nfevalmax - nfeval <= N_eps_relaxation*OneOtherLoop && ...
                        all(f_constraint_record > 0) ) || ...
                        flag_relaxation == 1 ) && flag_accept_relaxation == 1
                    
                    %                             any(f_inner_constraint > 0 + Delta_f_inner) && ...
                    %                             nfeval/nfevalmax > 0.8                     && ...
                    %                             all(Archive_matrix_constraint(:, problem_minmax.dim_d + 2*problem_minmax.dim_u + 3)~=0)
                    
                    flag_relaxation = 1;
                    
                    
                    f_constraint_record = [f_constraint_record; f_inner_constraint];
                    u_record_maxC{obj} = [u_record_maxC{obj}; umax_constraint];
                    
                    u_record_maxC{obj}  = u_record_maxC{obj} (f_constraint_record >= Delta_f_inner_min + (Delta_f_inner_max - Delta_f_inner_min)/N_eps_relaxation*N_reduce_record_constraint, :);
                    f_constraint_record = f_constraint_record(f_constraint_record >= Delta_f_inner_min + (Delta_f_inner_max - Delta_f_inner_min)/N_eps_relaxation*N_reduce_record_constraint, :);
                    
                    N_reduce_record_constraint = N_reduce_record_constraint + 1;
                end
                
            end
            %----------------------------------------------------------
            
            N_u_mpaidea = N_u_mpaidea + 1;
            
            
            % Archive best between umax and umax2 [in case IDEA failed]
            if any(sign_inner*f_inner > sign_inner*f_outer(i,obj))
                
                
                % if the maximisation in u give a higher (worse) F
                % sign_inner = 1 for minmax
                % f_inner = max F(dmin,u)
                % f_outer = min F(d,Au)
                
                if lsflag_inner
                    u_cell_aux = cell(1,n_obj);
                    u_cell_aux{obj} = umax;
                    [fmax2,umax2,nfeval2] = u_validation(problem_max_u, dmin(i,:), u_cell_aux, lsflag_inner, obj);
                    nfeval = nfeval + nfeval2;
                    
                    u_record{obj} = [u_record{obj}; umax2{obj}];
                    f_record_aux(i,obj) = fmax2(1,obj);
                else
                    u_record{obj} = [u_record{obj}; umax];
                    
                    
                    
                    dmin_true = dmin(i,:).*(problem_minmax.ub_d'-problem_minmax.lb_d') + (problem_minmax.lb_d');
                    umax_true = map_affine(umax, map_info);
                    
                    
                    
                    % save results in the Archive Table
                    [ARCHIVE, ARCHIVE_normalised] = archive_archive_normalised(ARCHIVE, ARCHIVE_normalised, problem_max_u, map_info, ...
                        umax_constraint, dmin_true, umax_true, f_inner, violation_u_kk, ...
                        f_inner_constraint, nfeval, dmin, umax, obj);
                    
                    
                    
                    
                    clearvars f_record_aux
                    %                         f_record_aux(i,obj) = f_inner;
                    f_record_aux = f_inner;
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
                    
                else
                    for k=1:size(u_record_tot{obj},1)
                        
                        d = dmin.*(problem_minmax.ub_d'-problem_minmax.lb_d') + (problem_minmax.lb_d');
                        map_info = problem_max_u.par_objfun.map_u_info{obj};
                        u = map_affine(u_record_tot{obj}(k,:),map_info);
                        
                        
                    if ~problem_max_u.par_objfun.obj_constr
                        f_outer_compare = problem_max_u.par_objfun.objfun{obj}(d, u, problem_max_u.par_objfun.problem_par_objfun{obj});
                    else
                        [out_obj_constr_compare] = problem_minmax.objfun{obj}(d, u, problem_max_u.par_objfun.problem_par_objfun{obj});
                        f_outer_compare = out_obj_constr_compare.f;
                    end
                        
                        
                        if f_outer(i) == f_outer_compare
                            
                            
                            % save in the Archive Table
                            [ARCHIVE, ARCHIVE_normalised] = archive_archive_normalised(ARCHIVE, ARCHIVE_normalised, problem_max_u, map_info, ...
                                umax_constraint, d, u, f_outer(i), violation_d_kk, ...
                                f_inner_constraint, nfeval, dmin, u_record_tot{obj}(k,:), obj);
                            
                            
                        end
                        nfeval = nfeval + 1;
                    end
                end
            end
            
            if algo_inner.save_archive
                
                fprintf(fileID_archive,str_archive,ARCHIVE{obj}(end,:));
                fprintf(fileID_archive_inner, str_archive_inner, ARCHIVE_INNER{obj}(end,:));
                fprintf(fileID_archive_outer, str_archive_outer, ARCHIVE_OUTER{obj}(end,:));
            end
        end
    end
    
    % Archive dmin
    d_record = [d_record; dmin];
    f_record = [f_record; f_record_aux];
    
    % Remove solutions archived more than once (if any)
    [~, idx] = unique(round(1e8*d_record),'rows');
    d_record = d_record(idx,:);
    f_record = f_record(idx,:);
    for obj = 1:n_obj
        [~, idx] = unique(round(1e8*u_record{obj}),'rows');
        u_record{obj} = u_record{obj}(idx,:);
        
        [~, idx] = unique(round(1e8*u_record_maxC{obj}),'rows');
        u_record_maxC{obj} = u_record_maxC{obj}(idx,:);
        
        f_constraint_record = f_constraint_record(idx);
    end
    
    
    % take into accout the function evaluation for the crosscheck
    size_u_record      = sum(cellfun('size',u_record,1));
    size_u_record_maxC = sum(cellfun('size',u_record_maxC,1));
    size_u_record_tot  = sum(cellfun('size',u_record,1)) + sum(cellfun('size',u_record_maxC,1));
    nfeval_loop = nfevalmax - size_u_record*size(d_record,1)*(~lsflag_validation + lsflag_validation*20*problem_minmax.dim_u);
    % NOTE: 50*dim_u is hardcoded nfevalmax in local search but is too conservative and overestimates nfeval_val hence 20.
    
    
    stop = nfeval >= nfeval_loop;
    save(strcat('nu400_20000inout_uncnstr_1obj_nfeval',num2str(nfeval)))
end

if algo_inner.save_archive
    
    fclose(fileID_archive);
    fclose(fileID_archive_inner);
    fclose(fileID_archive_outer);
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

%----------------------------------------------------------
% CONSTRAINTS
%
% no constraint

if isempty(problem_max_u.constraint) && ~problem_max_u.par_objfun.obj_constr
    
    while ~stop
        sel = dominance(f_record,0) == 0;       % find non-dominated
        if(n_obj>1)
            sel = sel';
        end
        
        tocheck = sel & ~checked;               % select only those that have not been checked
        stop = all(~tocheck);                   % stop if you have nothing to check
        % check those that have been selected
        d_tocheck = d_record(tocheck,:);
        
        if stop
            break
        end
        
        [f_val_aux, u_val_record_aux, nfeval_aux] = u_validation(problem_max_u, d_tocheck, u_record, lsflag_validation, 1:n_obj);
        
        
        nfeval_val = nfeval_val + nfeval_aux;
        
        
        % update f_record and u_val_record
        if isstruct(f_val_aux)
            f_record(tocheck,:) = f_val_aux.f;
        else
            f_record(tocheck,:) = f_val_aux;
        end
        
        for obj = 1:n_obj
            u_val_record{obj}(tocheck,:) = u_val_record_aux{obj};
        end
        checked = checked | tocheck;           % update who has been checked
    end
    fval = f_record(sel,:);
    nfeval = nfeval + nfeval_val;
    
    
    frontsize = size(fval,1);
    
    d = d_record(sel,:).*repmat(problem_minmax.ub_d'-problem_minmax.lb_d',[frontsize,1]) + repmat(problem_minmax.lb_d',[frontsize,1]);
    u = cell(1,n_obj);
    for obj = 1:n_obj
        u_val_record{obj} = u_val_record{obj}(sel,:);
        map_info = problem_max_u.par_objfun.map_u_info{obj};
        for i = 1:frontsize
            u{obj}(i,:) = map_affine(u_val_record{obj}(i,:),map_info); %this can be easily vectorized
        end
    end
    
    output.u = u;
    output.nfeval = nfeval;
    output.ARCHIVE = ARCHIVE;
    output.ARCHIVE_C = ARCHIVE_C;
    output.ARCHIVE_INNER = ARCHIVE_INNER;
    output.ARCHIVE_OUTER = ARCHIVE_OUTER;
    exitflag = 0;
    
    
    
else
    
    
    
    % -----------------------------------------------------------------
    % cross-check
    %         ARCHIVE_checked = ARCHIVE;
    %         U_check = [cell2mat(ARCHIVE{obj}(2:end,2)); cell2mat(ARCHIVE_C{obj}(2:end,2))];
    %         for N_d_archive = 2:size(ARCHIVE{1, 1},1)
    %             for N_u_archive = 1:size(ARCHIVE{1, 1},1)-1 + size(ARCHIVE_C{1, 1},1)-1
    %                     d_archive = ARCHIVE{1, 1}{N_d_archive, 1};
    %                     f_checked = problem_minmax.objfun{1}(d_archive, U_check(N_u_archive, :), []);
    %                     c_checked = problem_minmax.constraints{1}(d_archive, U_check(N_u_archive, :), []);
    %                     if f_checked > ARCHIVE{1, 1}{N_d_archive, 3} && c_checked <= ARCHIVE{1, 1}{N_d_archive, 4}
    %                         ARCHIVE_checked{1,1}{N_d_archive, 2} = U_check(N_u_archive, :);
    %                         ARCHIVE_checked{1,1}{N_d_archive, 3} = f_checked;
    %                     end
    %             end
    %         end
    
    
    f_archive = cell2mat(ARCHIVE{obj}(2:end,4));
    c_archive = cell2mat(ARCHIVE{obj}(2:end,5));
    
    ARCHIVE_crosscheck = ARCHIVE_normalised;
    ARCHIVE_crosscheck{1, 1}(:,1) = ARCHIVE{1, 1}(:,1);
    [ARCHIVE_crosscheck, nfeval_aux] = u_validation_constraint(problem_max_u, c_archive, lsflag_validation, 1:n_obj, ARCHIVE_normalised, ARCHIVE_crosscheck);
    nfeval = nfeval + nfeval_aux;
    % -----------------------------------------------------------------
    
    
    %         pop_inner_con = find(ARCHIVE{obj}(:, problem_minmax.dim_d + problem_minmax.dim_u + 2) < EPS_constraint);
    %         pop_max_con   = find(ARCHIVE{obj}(:, problem_minmax.dim_d + problem_minmax.dim_u + 3) < EPS_constraint);
    Archive_matrix = cell2mat(ARCHIVE_crosscheck{obj}(2:end,:));
    pop_inner_con = find(Archive_matrix(:, problem_minmax.dim_d + 2*problem_minmax.dim_u + 2) <= EPS_constraint);
    pop_max_con   = find(Archive_matrix(:, problem_minmax.dim_d + 2*problem_minmax.dim_u + 3) <= EPS_constraint);
    
    %         pop_inner_con = find(ARCHIVE{obj}{2:end, 4} < EPS_constraint);
    %         pop_max_con   = find(ARCHIVE{obj}{2:end, 5} < EPS_constraint);
    
    feasible = intersect(pop_inner_con, pop_max_con);
    if ~isempty(feasible)
        ARCHIVE_selected{obj} = Archive_matrix(feasible, :);
        
    elseif ~isempty(pop_inner_con)
        ARCHIVE_selected{obj} = Archive_matrix(pop_inner_con, :);
        
    else
        NOTpop_inner_con = find(Archive_matrix(:, problem_minmax.dim_d + 2*problem_minmax.dim_u + 2) ~= 0);
        NOTpop_max_con   = find(Archive_matrix(:, problem_minmax.dim_d + 2*problem_minmax.dim_u + 3) ~= 0);
        NOTfeasible = intersect(NOTpop_inner_con, NOTpop_max_con);
        ARCHIVE_selected{obj} = Archive_matrix(NOTfeasible, :);
        
    end
    
    [f_minmax, pos_minmax] = min(ARCHIVE_selected{obj}(:, problem_minmax.dim_d + 2*problem_minmax.dim_u + 1));
    d = ARCHIVE_selected{obj}(pos_minmax, 1:problem_minmax.dim_d);
    u = ARCHIVE_selected{obj}(pos_minmax, problem_minmax.dim_d+1:problem_minmax.dim_d + problem_minmax.dim_u);
    %         d = d_scaled.*(problem_minmax.ub_d'-problem_minmax.lb_d') + (problem_minmax.lb_d');
    %         map_info = problem_max_u.par_objfun.map_u_info{obj}; output.u{obj}   = map_affine(u_scaled, map_info);
    output.u{obj} = u;
    output.nfeval = nfeval;
    output.ARCHIVE = ARCHIVE;
    output.ARCHIVE_crosscheck = ARCHIVE_crosscheck;
    output.ARCHIVE_C = ARCHIVE_C;
    output.ARCHIVE_INNER = ARCHIVE_INNER;
    output.ARCHIVE_OUTER = ARCHIVE_OUTER;
    exitflag = 0;
    fval = f_minmax;
    
    
    
    %         for d_index = 1:size(d_record)
    %
    %             d_tocheck = d_record(d_index,:);
    %
    %             [f_val_aux, u_val_record_aux, nfeval_aux, ~, violation] = u_validation_constraints(problem_max_u, d_tocheck, u_record, lsflag_validation, 1:n_obj);
    %
    %             nfeval_val = nfeval_val + nfeval_aux;
    %
    %
    %             % update f_record and u_val_record
    %             f_record(d_index,:) = f_val_aux;
    %             f_violation(d_index,:) = violation;
    %
    %             for obj = 1:n_obj
    %                 u_val_record{obj}(d_index,:) = u_val_record_aux{obj};
    %             end
    %
    %         end
    %
    %         sel = dominance(f_record,0) == 0;       % find non-dominated
    %         if f_violation(sel) ~= 0
    %
    %             if any(f_violation==0)
    %                 [m_f, ~] = min(f_record(f_violation==0));
    %                 sel = find(f_record == m_f);
    %             else
    %                 [m_c, ~] = min(f_violation(f_violation>0));
    %                 sel = find(f_violation == m_c);
    %             end
    %         end
    
    
    
end







end

% end