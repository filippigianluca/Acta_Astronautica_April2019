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

%% ------------------------------------------------------------------------
% Reset random numbers generator
%--------------------------------------------------------------------------
seed = 1;
s = RandStream('mt19937ar','Seed', seed);
RandStream.setGlobalStream(s);




%% ------------------------------------------------------------------------
% Define Problem
%--------------------------------------------------------------------------

[problem, algo_inner, algo_outer, algo_minmax, algo_decomposition] = def_ebro_IAC_zeno_SA();
% [problem, algo_inner, algo_outer, algo_minmax, algo_decomposition] = def_ebro_IAC();




%% ------------------------------------------------------------------------
% RUN MINMAX and/or MINMIN PROBLEM
%--------------------------------------------------------------------------

% for i = 5000:5000:40000
% algo_inner.par.nFeValMax = i;

[minmin, minmax, LIST, LIST_EXACT] = ebro(problem, algo_inner, algo_outer, algo_minmax, algo_decomposition);

% save(strcat('zeno_sensitivity_max_Dfix',num2str(i)))
% end
