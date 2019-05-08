% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [EXACT_FE, LIST_EXACT] = evaluate_Belief_Plausibility_exact(in, minmax, minmin, n_obj, algo_decomposition)

%--------------------------------------------------------------------------
% FREEZE U AND D FROM INPUT
%--------------------------------------------------------------------------
problem_inner = build_metaproblem_ideaminmax_s_inner(in);

if in.flag_output.exact_Belief
    problem_inner.par_objfun.d_belief = minmax.d;
    problem_inner.par_objfun.u_belief = minmax.u;
end

if in.flag_output.exact_Plausibility
    problem_inner.par_objfun.d_plausibility = minmin.d;
    problem_inner.par_objfun.u_plausibility = minmin.u;
end





%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
% problem_inner.par_objfun.d = dmin;
problem_inner.par_objfun.objective = n_obj;
% problem_inner.par_objfun.umax = u_max_tot;
% problem_inner.objfun = @mask_objfun_max_decomposition;
problem_inner.par_objfun.objfun = in.objfun;
problem_inner.constraint = in.constraint{n_obj};
problem_inner.par_objfun.problem_par_objfun{n_obj} = in.par_objfun{n_obj};



% structure decomposition
num_FE_TOT = number_Focal_Element_TEST(in);
EXACT_FE = cell(1,num_FE_TOT);




problem_inner.par_objfun.vars_to_opt = [1:problem_inner.dim];              % all the components of vector u
problem_inner.dim = length(problem_inner.par_objfun.vars_to_opt);

% init
problem_max = in;

problem_max.dim_u = problem_inner.dim;


for index =1:problem_inner.dim
    lb_u(index,1)      = {in.lb_u{n_obj}{index}};
    ub_u(index,1)      = {in.ub_u{n_obj}{index}};
    bpa(index,1)       = {in.bpa{n_obj}{index}};
    lengt_bpa(1, index) = length(bpa{index});
end
problem_max.lb_u = lb_u;
problem_max.ub_u = ub_u;
% BPA
problem_max.bpa = bpa;




% initialize algorithm
problem_inner.n_obj = n_obj;
problem_inner.dim_d = in.dim_d;
problem_inner.dim_u = sum(in.dim_u);


problem_inner.bpa = problem_max.bpa;

in.dim_u = problem_inner.dim;

pos_FE = '[position_FE(1)';
for iii=2:length(lb_u)
    pos_FE = strcat(pos_FE, ',position_FE(',num2str(iii),')'); 
end
pos_FE = strcat(pos_FE, ']'); 

LIST_EXACT.FE        = [];
LIST_EXACT.FE_Belief = {'N FE', 'F_Bel', 'u_Bel', 'bpa_Bel', 'LB', 'UB'};



for FE_i = 1:num_FE_TOT
    problem_inner.lb = zeros(1,problem_inner.dim);
    problem_inner.ub = zeros(1,problem_inner.dim);
    
    

%     position_FE = position(ii, problem_inner, problem_max);   
%     eval(strcat(pos_FE, ' = ind2sub(3*ones(1,',num2str(problem_max.dim_u),'),', num2str(FE_i), ');'));   
    eval(strcat(pos_FE, ' = ind2sub([',num2str(lengt_bpa),'],', num2str(FE_i), ');'));


    for iii = 1 : problem_inner.dim
        
        problem_inner.lb(iii) = problem_max.lb_u{iii,1}(position_FE(iii));
        problem_inner.ub(iii) = problem_max.ub_u{iii,1}(position_FE(iii));
        
    end

%     A= 1;
%     for L = 0.1:0.1:1
%         for M = 0.1:0.1:1
%         algo_decomposition.par.F = [];
%         algo_decomposition.par.CR = [];    
    
    [EXACT_FE] = max_func_decomposizion_TEST(in, FE_i, position_FE, EXACT_FE, algo_decomposition, problem_inner);
%     EXACT_FE{1, 1}.N= A;
%     EXACT_FE{1, 1}.L= L;
%     EXACT_FE{1, 1}.CR= M;
%     save(strcat('IAC_max_exact_1000feval',num2str(A),'10LRmax'))
%     A=A+1;
%         end
%     end
    
    
    if in.flag_output.exact_Belief
        
            LIST_EXACT.FE_Belief = [LIST_EXACT.FE_Belief;...
                                 FE_i...
                                 EXACT_FE{FE_i}.upper_f...
                                 {EXACT_FE{FE_i}.upper_u}...
                                 EXACT_FE{FE_i}.bpa...
                                 {problem_inner.lb}...
                                 {problem_inner.ub}...
                            ];
    end
    
    if in.flag_output.exact_Plausibility
            LIST_EXACT.FE_Plausibility = [LIST_EXACT.FE;...
                                         FE_i...
                                         EXACT_FE{FE_i}.downer_f...
                                         EXACT_FE{FE_i}.downer_u...
                                         EXACT_FE{FE_i}.bpa...
        %                      {problem_inner.lb}...
        %                      {problem_inner.ub}...
                            ];
    end
    
    
    %             savefolder = strcat('RESULTS_spacecraft_1_16000/exact/');
    %             mkdir(savefolder);
    %             save(strcat(savefolder,'spacecraft_exact_1_16000'));
    
    
end


global num_maximization_exact;
num_maximization_exact = num_FE_TOT;


%% plot Bl

for k = 1 : length(EXACT_FE)     % vector of all the f-values of the FEs
    if in.flag_output.exact_Belief
        f_belief(k) = EXACT_FE{k}.upper_f;
    end
    
    if in.flag_output.exact_Plausibility
        f_plausibility(k) = EXACT_FE{k}.downer_f;
    end
    
    bpa_end(k) = EXACT_FE{k}.bpa;
end


if in.flag_output.exact_Belief
    
    [f_sorted_Bel, position_sorted_Bel] = sort(f_belief);
    
    LIST_EXACT.F_Bel = [f_sorted_Bel(1), f_sorted_Bel];
    LIST_EXACT.Bel   = [0, cumsum(bpa_end(position_sorted_Bel))];
end


% Plausibility
if in.flag_output.exact_Plausibility
    
    [f_sorted_Pl, position_sorted_Pl] = sort(f_plausibility);
    
    LIST_EXACT.F_Pl = [f_sorted_Pl(1), f_sorted_Pl];
    LIST_EXACT.Pl   = [0, cumsum(bpa_end(position_sorted_Pl))];
end




end