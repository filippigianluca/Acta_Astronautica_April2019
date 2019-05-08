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


warning('off','all')
%warning('on','all')
%--------------------------------------------------------------------------
% add folders to path
%--------------------------------------------------------------------------
CURRENTPATH=pwd;
if isunix
    idcs   = strfind(CURRENTPATH,'/');
else
    idcs   = strfind(CURRENTPATH,'\');
end

ALLdir = CURRENTPATH(1:idcs(end-2)-1);
IAC2018 = CURRENTPATH(1:idcs(end)-1);
if isunix
    IAC2018_END = strcat(IAC2018,'/FORMULATION_END');
else
    IAC2018_END = strcat(IAC2018,'\FORMULATION_END');
end
rmpath(genpath(ALLdir));
addpath(genpath(IAC2018_END));


% % % ALLdir1 = CURRENTPATH(1:idcs(end-2)-1);
% % % IAC2018 = CURRENTPATH(1:idcs(end)-1);
% % % Optimisation_smarto2c_dir = strcat(ALLdir1,'\SMART_O2C\smart-o2c\Optimisation');
% % % Problem_smarto2c_dir = strcat(ALLdir1,'\SMART_O2C\smart-o2c\Problems');
% % % IAC_formulation_1_dir     = strcat(ALLdir1,'\UTOPIAE_RELIABLES\IAC_2018\FORMULATION_1');
% % % IAC_formulation_2_dir     = strcat(ALLdir1,'\UTOPIAE_RELIABLES\IAC_2018\FORMULATION_2');
% % % 
% % % addpath(genpath(IAC2018));
% % % addpath(genpath(Optimisation_smarto2c_dir));
% % % 
% % % rmpath(genpath(Problem_smarto2c_dir));
% % % rmpath(genpath(IAC_formulation_1_dir));
% % % rmpath(genpath(IAC_formulation_2_dir));


%--------------------------------------------------------------------------
% Reset random numbers generator
%--------------------------------------------------------------------------
seed = 1;
s = RandStream('mt19937ar','Seed', seed);
RandStream.setGlobalStream(s);


%--------------------------------------------------------------------------
% Define Problem
%--------------------------------------------------------------------------
[problem] = def_IAC_2018();
[algo_inner, algo_outer, algo_minmax, algo_decomposition] = def_optimiser_IAC_2018_so(problem);


% [problem] = init_tc_so_constr_IAC_2018();
% [algo_minmax, algo_outer, algo_inner, algo_decomposition] = init_algo_so_ebro_IAC_2018(problem);



global Mass
global Mpayload 
global Mobdh 
global Maocs 
global Mttc 
global Mpower;
global d_ 
global u_

d_ = [];
u_ = [];
Mass = [];
Mpayload = [];
Mobdh = [];
Maocs = [];
Mttc = [];
Mpower = [];

[minmin, minmax, LIST, LIST_EXACT] = ebro(problem, algo_inner, algo_outer, algo_minmax, algo_decomposition);

save('f_minmax_10e6_2e3_2e3')
