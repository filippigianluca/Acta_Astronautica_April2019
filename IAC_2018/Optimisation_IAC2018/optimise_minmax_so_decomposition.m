% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [best_case, worst_case, LIST, LIST_EXACT] = optimise_minmax_so_decomposition(problem, algo_inner, algo_outer, algo_minmax, algo_decomposition)
% =========================================================================
% EBRO (Evidence-Based Robust Optimisation):
% - UQ using epistemic uncertainty;
% - Evidence theory with Belief and Plausibility curves;
% - Decomposition approach to treat with complex systems: uncertain
%   variables that couple two or more sub-systems and uncertain variables
%   that influence only one subsystems are studied separatedly.
% - The Decomposition approach assure conservative curves at  polynomial
%   cost under some assumptions (see the reference papers).
%
% Reference:
% M. Vasile, G. Filippi, C. O. Absil, A. Riccardi, "Fast Belief Estimation
% in Evidence Network Models", EUROGEN 2017, September 13-15,Madrid, Spain.
% G. Filippi, M. Vasile, M. Marchi, P. Vercesi, "Evidence-Based Robust
% Optimisation of Space Systems with Evidence Network Models", IEEE World
% Congress on Computational Intelligence, 8-13 July 2018, Rio de Janeiro,
% Brazil.
% =========================================================================
% =========================================================================
%
% INPUT
%
% * problem: structure containing the informations of the problem:
%          * problem.output: choose to reconstruct Belief (0),  Plausibility
%            (1) or both Belief and Plausibility (2);
%          * problem.input: run minmax and/or minmin (0), load the design
%            vector d (1) or load the design vector d and the uncertain
%            vector u (2);
%          * problem.exact_curves: reconstruct Belief and/or Plausibility
%            curves with the decomposition approach but do not evaluate the
%            exact curve(s) (0), do bot decomposition and exact curve (1)
%            or do only exact curve(s) (2);
%          * problem.num_functions: number of sub-functions in which the
%            system can be decomposed;
%          * problem.num_samples: number of samples in the partial curves;
%            the bigger the more precise the decomposition curve;
%          * problem.dim_u: dimention of the uncertain space [u1, u2, u12];
%          * problem.lb_u{i}: lower boundaries of the uncertain variables
%            for the i-th objective function;
%          * problem.ub_u{i}: upper boundaries of the uncertain variables
%            for the i-th objective function;
%          * problem.bpa: basic probability assignment of all the intervals
%            of each epistemic variable;
%          * problem.dim_d: dimention of the design space;
%          * problem.lb_d: lower boundaries of the uncertain variables
%          * problem.ub_d: upper boundaries of the uncertain variables
%          * problem.fix: structure of the fixed parameters (if any);
%          * problem.n_obj: number of objective function;
%          * problem.objfun = {@...}: objective function;
%          * problem.constraints = {@...}: constraint function;
%
% * algo_minmax: structure containing MP-AIDEA specific input information
%                (if empty default values are used) for the meta-algorithm.
%              * algo_minmax.par_minmax.maxnfeval: max number of function
%                evaluation for the minmax problem;
%
% * algo_outer: structure containing MP-AIDEA specific input information
%                (if empty default values are used) for the minimisation loop.
%              * algo_outer.par_mpaidea.nFeValMax: max number of function
%                evaluation for the outer loop (minimisation over d);
%              * algo_outer.par_mpaidea.n_populations;
%              * algo_outer.par_mpaidea.n_agents:
%              * algo_outer.par_mpaidea.max_LR: maximum number of local
%                restart before global restart (for cases when only one
%                population is considered and no adaptation of delta_local and
%                local/global restart is performed; if par_mpaidea.max_LR = []
%                 then adaptation is performed);
%
% * algo_inner: structure containing MP-AIDEA specific input information
%                (if empty default values are used) for the minimisation loop.
%             * algo_inner.par_mpaidea.nFeValMax;
%             * algo_inner.par_mpaidea.n_populations;
%             * algo_inner.par_mpaidea.n_agents;
%             * algo_inner.par_mpaidea.max_LR;
%
%
%
%% OUTPUT
%
% * minmax: structure with the worst case solution:
%         * minmax.d;
%         * minmax.u;
%         * minmax.f;
%
% * minmin: structure with the best case solution:
%         * minmin.d;
%         * minmin.u;
%         * minmin.f;
%
% * LIST: structur with the information of the Focal Elements considered in
%         the decomposition approach;
%
% * LIST EXACT: structur with the information of all the Focal Elements


best_case  = [];
worst_case = [];
best_case_ebro  = [];
worst_case_ebro = [];
LIST       = [];
LIST_EXACT = [];
problem_ebro = problem;
% n_obj_tot = problem.n_obj;



%% ------------------------------------------------------------------------
% DO min-max and/or min-min
% -------------------------------------------------------------------------

if problem.flag_output.minmin || problem.flag_output.minmax
    
    [ best_case, worst_case ] = optimise_minmax_so(problem, algo_outer, algo_inner, algo_minmax.par_minmax);
    
end




%% ------------------------------------------------------------------------
% load d and/or u and/or f for the minimisation problem
% -------------------------------------------------------------------------

if problem.flag_input.dmin_load
    best_case.d = problem.input.d_min;
end

if problem.flag_input.umin_load
    best_case.u = problem.input.u_min;
end

if problem.flag_input.dmin_load && problem.flag_input.umin_load && problem.flag_input.fmin_load
    best_case.f = problem.input.F_min;
elseif problem.flag_input.dmin_load && problem.flag_input.umin_load && ~problem.flag_input.fmin_load
    best_case.f = problem.objfun(worst_case.d, worst_case.u{1}, problem.fix);
end





%% ------------------------------------------------------------------------
% load d and/or u and/or f for the maximisation problem
% -------------------------------------------------------------------------

if problem.flag_input.dmax_load
    worst_case.d = problem.input.d_max;
end

for n_obj = 1:problem_ebro.n_obj
    
    if problem.flag_input.umax_load
        worst_case.u{n_obj} = problem.input.u_max{n_obj};
    end
    
    if problem.flag_input.dmax_load && problem.flag_input.umax_load && problem.flag_input.fmax_load
        worst_case.f{n_obj} = problem.input.F_max{n_obj};
    elseif problem.flag_input.dmax_load && problem.flag_input.umax_load && ~problem.flag_input.fmax_load
        worst_case.f{n_obj} = problem.objfun{n_obj}(worst_case.d, worst_case.u{n_obj}, problem.par_objfun{n_obj} );
    end
    
end





%% ------------------------------------------------------------------------
% DO maximisation and/or minimisation for fixed design vector(s)
% -------------------------------------------------------------------------

% to extend with the MO
if problem.flag_output.u_run
    
    if problem.flag_output.Belief
        problem.sign_inner = 1;
        [worst_case] = evaluate_max(problem, worst_case.d, algo_inner);
        
    end
    
    if problem.flag_output.Plausibility
        problem.sign_inner = -1;
        [best_case] = evaluate_min(problem, best_case.d, algo_inner);
        
    end
    
end



%% ------------------------------------------------------------------------
% DO a minimisation over d for a fixed uncertain vector
% -------------------------------------------------------------------------
if problem.flag_output.dmin_run
    
    minmax  = [];
    minmin  = [];
    maximum = [];
    minimum = [];
    
    problem.par_objfun{1, 1}.fix = problem.par_objfun{1}.fix;
    problem.par_objfun{1, 1}.fix.nominal_u = problem.input.vector_u_nominal;
    problem.par_objfun{1, 1}.objfun = problem.objfun;
    
    problem.dim = length(problem.lb_d);
    problem.lb = problem.lb_d';
    problem.ub = problem.ub_d';
    
    problem.objfun = @mask_minimise_d_fixed_u;
    problem.constraint = problem.constraint{1};
    
    
%     problem.fitnessfcn.obj        = @mask_minimise_d_fixed_u;
%     problem.fitnessfcn.constr     = problem.constraint{1};
%     problem.fitnessfcn.obj_constr = 0;
%     problem.fitnessfcn.weighted   = 0;
%     problem.fitnessfcn.ceq_eps    = 1e-6;
%     problem.fitnessfcn.w_ceq      = 100;
%     problem.fitnessfcn.w_c        = 100;
    
    
    [minimum] = algo_outer.optimise(problem, algo_outer.par);
    min_F = problem.par_objfun{1, 1}.objfun{1, 1}  (minimum, problem.input.vector_u_nominal, problem.par_objfun{1});
    maximum_new = [];
    
end





%% ------------------------------------------------------------------------
% DO decomposition approach for Belief and/or Plausibility reconstruction,
% or evaluate exact Belief and/or exact Plausibility
% -------------------------------------------------------------------------


if problem.flag_output.Belief || problem.flag_output.Plausibility || problem.flag_output.exact_Belief || problem.flag_output.exact_Plausibility
    
    for n_obj = 1:problem_ebro.n_obj
        
        if problem.flag_output.Belief || problem.flag_output.exact_Belief
            worst_case_ebro.d = worst_case.d;
            worst_case_ebro.u = worst_case.u{n_obj};
            worst_case_ebro.f = worst_case.f{n_obj};
        end
        if problem.flag_output.Plausibility || problem.flag_output.exact_Plausibility
            best_case_ebro.d = best_case.d;
            best_case_ebro.u = best_case.u{n_obj};
            best_case_ebro.f = best_case.f{n_obj};
        end
        
        [best_case, worst_case, LIST, LIST_EXACT] = ebro_decomposition(problem_ebro, algo_decomposition, worst_case_ebro, best_case_ebro);
        
    end
end





return
