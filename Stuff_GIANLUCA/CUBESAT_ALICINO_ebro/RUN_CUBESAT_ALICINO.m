% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of run of a min-max problem using MP_AIDEA in a single objective
% optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all; clc

%--------------------------------------------------------------------------
% Reset random numbers generator
%--------------------------------------------------------------------------
seed = 1;
s = RandStream('mt19937ar','Seed', seed);
RandStream.setGlobalStream(s);


%--------------------------------------------------------------------------
% Define Problem
%--------------------------------------------------------------------------

[problem, algo_inner, algo_outer, algo_minmax, algo_decomposition] = def_CUBESAT_ALICINO();





%% ------------------------------------------------------------------------
% RUN MINMAX and/or MINMIN PROBLEM
%--------------------------------------------------------------------------
% 
% for k = 1:10
%     MAX = [];
%     problem.input_d_max  = problem.lb_d + rand(9,1).*(problem.ub_d - problem.lb_d);
% for i = [1000:1000:10000]
%     algo_inner.par.nFeValMax = i;
    
    
% % 
% %     for L =0.1:0.1:1
% %         for M = 0.1:0.1:1
%             
%             algo_inner.par.F = [];                           % F
%             algo_inner.par.CR = [];                          % CR
    
    
    [minmin, minmax, LIST, LIST_EXACT] = ebro(problem, algo_inner, algo_outer, algo_minmax, algo_decomposition);

%     save('exact_Belief_CUBESAT_alicino_paper_29luglio_5000evalFE')

% MAX = [ MAX minmax.f ];    
% save(strcat('CUBESAT_MAX_ALICINO_29LUGLIO',num2str(k)))
% save(strcat('IAC_max_exact',num2str(i),'_feval_20LR_CRandF_adaptive_025deltalocal'))
% %     
% %         end
% %     end
% end
% end
