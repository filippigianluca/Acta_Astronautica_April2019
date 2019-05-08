% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [minmin, minmax, LIST, LIST_EXACT] = ebro_decomposition(problem, algo_decomposition, minmax_input, minmin_input)

disp(strcat('Running Decomposition ...'));

n_points_tot = 1;
n_obj = problem.n_obj;

minmax     = [];
minmin     = [];
LIST       = [];
LIST_EXACT = [];



for n_point = 1:n_points_tot                  % all the design points(DP); if single objective there is only one DP.
    
    
    
    if isfield(minmax_input,'d') && ~isempty(minmax_input.d)
        minmax.d = minmax_input.d(n_point,:);
        minmax.u = minmax_input.u(n_point,:);
        minmax.f = minmax_input.f(n_point, n_obj);
        
        minmax.FE = find_FE(n_obj, minmax.u, problem);
    end
    if isfield(minmin_input,'d') && ~isempty(minmin_input.d)
        minmin.d = minmin_input.d(n_point,:);
        minmin.u = minmin_input.u(n_point,:);
        minmin.f = minmin_input.f(n_point, n_obj);
        
        minmin.FE = find_FE(n_obj, minmin.u, problem);
    end
    
    
    
     
    problem.sign_inner = 1;
    
    if problem.flag_output.Belief || problem.flag_output.Plausibility       
        %% NETWORK ANALYSIS
        problem = network_analysis(problem, minmax, minmin);
        
        
        %% PARTIAL CURVES (coupled variables)
        [decomposition, Partial_curve] = evaluate_Belief_Plausibility_coupled_vectors(problem, minmax, minmin, n_obj, n_point, algo_decomposition);
        
        
        %% SAMPLE
        [Sample] = Sampling_Belief_Plausibility_coupled_u(Partial_curve, decomposition, problem, minmax, minmin);
        
        
        %% RECONSTRUCTION
        [decomposition_end, Plot_decomposition, num_sample_tot, LIST] = reconstruction_Belief_Plausibility(problem, minmax, minmin, Sample, n_obj, n_point, algo_decomposition, Partial_curve);
        
        
    end
    
    
    
    
    %% EXACT CURVES (Belief and/or Plausibility)
    if problem.flag_output.exact_Belief || problem.flag_output.exact_Plausibility
        
        [EXACT_FE, LIST_EXACT] = evaluate_Belief_Plausibility_exact(problem, minmax, minmin, n_obj, algo_decomposition);
        
    end
    
    
    if   problem.flag_output.plot
        
        plot_ebro_decomposition(problem, Partial_curve, LIST, LIST_EXACT);
    end
    
    
    
end


return