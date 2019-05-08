% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [d,fval,exitflag,output] = minmin_so(problem_minmax, algo_outer, algo_inner, par_minmax)


% Rename inputs
n_d = problem_minmax.dim_d;
n_u = problem_minmax.dim_u;
n_obj = problem_minmax.n_obj;


% 1 for minmax, -1 for minmin
sign_inner = -1;
problem_minmax.sign_inner = -1;


nfevalmax = par_minmax.maxnfeval;

N_u_mpaidea = 2;
N_txt_mpaidea_inner  = 1;
N_txt_mpaidea_constr = 1;
N_txt_mpaidea_outer  = 1;


EPS_constraint = 0;  % realmin('single');

% doesn't accept the relaxation
flag_accept_relaxation = 0;



if isempty(problem_minmax.constraint{1})
    
    %% --------------------------------------------------------------------
    % UN-CONSTRAINED MINMIN
    %----------------------------------------------------------------------
    
    % Define the Problem
    problem_minmin = build_metaproblem_minmin_so(problem_minmax);
    
    
    par = algo_inner.par;
    par.nFeValMax = nfevalmax;
    
    
    par.sign_inner=sign_inner;
    problem_minmin.par_objfun.sign = sign_inner;
    
%     problem_minmin.par_objfun.problem_par_objfun{n_obj}.fix = problem_minmax.fix;  % FIXED PARAMETERS
    
    
    % Run minimisation over all vector [d, u]
    [ dumin, fval_all_populations , exitflag , output_aux] = algo_inner.optimise(problem_minmin,par);
    
    
    
    [fval, min_pop] = min(fval_all_populations);
    dmin = dumin(min_pop, 1:n_d);
    umin = dumin(min_pop, n_d+1:end);
    d = dmin.*(problem_minmax.ub_d'-problem_minmax.lb_d') + problem_minmax.lb_d';
    u = cell(1,n_obj);
    obj=1;
    map_info = problem_minmin.par_objfun.map_u_info{obj};
    u{obj} = map_affine(umin,map_info);
    
    output.u = u;
    output.nfeval = output_aux.nfeval;
    
    
    
    
    
    
    
    
else
    %% --------------------------------------------------------------------
    % CONSTRAINED MINMIN
    %----------------------------------------------------------------------
    
    % Initialise Archives
    Archive_u_objfun = [];
    Archive_d_objfun = [];
    Archive_f_objfun = [];
    Archive_u_constr = [];
    
    % Define Problems
    problem_max_u_C = build_metaproblem_macsminmax_inner(problem_minmax);
    problem_minmin  = build_metaproblem_minmin_so(problem_minmax);
    
    nfeval = 0;
    
    stop = false;
    
    while ~stop
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % lower level: constrained minimisation of obj;
        % the onstraint function has to be satisfied in an archive of
        % uncertaintain vectors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        par = algo_inner.par;
        %         par.nFeValMax = nfevalmax;
        
        
        par.sign_inner=sign_inner;
        problem_minmin.par_objfun.sign = sign_inner;
        

        
        
        % FIXED PARAMETERS
        problem_minmin.par_objfun.problem_par_objfun{n_obj}.fix                          = problem_minmax.par_objfun{1, 1}.fix;   
        problem_minmin.par_objfun.problem_fix_d.par_objfun.problem_par_objfun{n_obj}.fix = problem_minmax.par_objfun{1, 1}.fix;   
        problem_minmin.par_objfun.problem_par_objfun{n_obj}.fix                          = problem_minmax.par_objfun{1, 1}.fix;   
        problem_minmin.par_objfun.problem_par_objfun{n_obj}.ub_u                         = problem_minmax.ub_u{1};
        problem_minmin.par_objfun.problem_par_objfun{n_obj}.lb_u                         = problem_minmax.lb_u{1};
        problem_minmin.par_objfun.problem_fix_d.par_objfun.lb_d                          = problem_minmax.lb_d';
        problem_minmin.par_objfun.problem_fix_d.par_objfun.ub_d                          = problem_minmax.ub_d';
        problem_minmin.par_objfun.problem_fix_d.par_objfun.constraint{1}                 = problem_minmax.constraint{1};
        problem_minmin.par_objfun.problem_fix_d.par_objfun.map_u_info{1}                 = problem_minmin.par_objfun.map_u_info{1};
        problem_minmin.par_objfun.objectives = 1;
        problem_minmin.par_objfun.u_record{1} = Archive_u_constr;
        problem_minmin.par_objfun.problem_fix_d.par_objfun.sign = -1;
        %%%
        
        [ dumin, fval_all_populations , exitflag , output_aux] = algo_outer.optimise(problem_minmin,par);
        
        nfeval = nfeval + output_aux.nfeval;
        
        [fval, min_pop] = min(fval_all_populations);
        dmin = dumin(min_pop, 1:n_d);
        umin = dumin(min_pop, n_d+1:end);
        
        
        dmin(dmin < 0) = 0;
        dmin(dmin > 1) = 1;
        umin(umin < 0) = 0;
        umin(umin > 1) = 1;
        
        
        d = dmin.*(problem_minmax.ub_d'-problem_minmax.lb_d') + problem_minmax.lb_d';
        u = cell(1,n_obj);
        obj=1;
        map_info = problem_minmin.par_objfun.map_u_info{obj};
        u{obj} = map_affine(umin,map_info);
        
        Archive_u_objfun = [Archive_u_objfun; u];
        Archive_d_objfun = [Archive_d_objfun; d];
        Archive_f_objfun = [Archive_f_objfun; fval];
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % upper level: maximisation of the constraint violation. The u
        % vector solution is then added to the Archive_u_objfun.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        problem_max_u_C.par_objfun.d = dmin;
        problem_max_u_C.par_objfun.sign = 1;
        problem_max_u_C.par_objfun.objective = 1;
        [ umax_constraint, f_inner_constraint , ~ , output_aux] = optimise_constraint(problem_max_u_C, algo_inner.par);
        
        nfeval = nfeval + output_aux.nfeval;
        
        umax_constraint(umax_constraint < 0) = 0;
        umax_constraint(umax_constraint > 1) = 1;
        
        if f_inner_constraint > 0
            Archive_u_constr = [Archive_u_constr; umax_constraint];
        end
        
        stop = (nfeval > par_minmax.maxnfeval );
        
    end


    [fval, pos_minmin_f] = min(Archive_f_objfun);

    d = Archive_d_objfun(pos_minmin_f, :);
    output.u = Archive_u_objfun(pos_minmin_f, :);
    output.nfeval = output_aux.nfeval;
    
    
    
end