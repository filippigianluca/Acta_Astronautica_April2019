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
clear all; close all; clc;



% add folders to path
CURRENTPATH=pwd;
if isunix
    idcs   = strfind(CURRENTPATH,'/');
else
    idcs   = strfind(CURRENTPATH,'\');
end
newdir = CURRENTPATH(1:idcs(end-1)-1); 
ALLdir = CURRENTPATH(1:idcs(end-1)-1); 

% remove everything
rmpath(genpath(ALLdir));



Gianlucadir=[newdir,'/Stuff_GIANLUCA'];
Dandadir=[newdir,'/Stuff_DANDA'];
newdir2=[newdir,'/Stuff_GIANLUCA/Optimisation_Multi_Objective'];
addpath(genpath(Gianlucadir));
addpath(genpath(Dandadir));
addpath(genpath(CURRENTPATH));
rmpath(genpath(newdir2));



%--------------------------------------------------------------------------
% Reset random numbers generator
%--------------------------------------------------------------------------
seed = 1;
s = RandStream('mt19937ar','Seed', seed);
RandStream.setGlobalStream(s);


%--------------------------------------------------------------------------
% Define Problem
%--------------------------------------------------------------------------

% functions from "A Surrogate-Assisted Evolutionary Algorithm for Minimax
% Optimization" by Aimin Zhou and Qingfu Zhang
%
% TC = [1,...,5] 


[problem] = def_SSTL();


[algo_inner, algo_outer, algo_minmax, algo_decomposition] = def_optimiser_SSTL(problem);



[minmin, minmax, LIST, LIST_EXACT] = ebro(problem, algo_inner, algo_outer, algo_minmax, algo_decomposition);


