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
    Stuff_Danda = strcat(ALLdir,'/UTOPIAE_RELIABLES/Stuff_DANDA');
else
    Stuff_Danda = strcat(ALLdir,'\UTOPIAE_RELIABLES\Stuff_DANDA');
end
rmpath(genpath(ALLdir));
addpath(genpath(IAC2018));
addpath(genpath(Stuff_Danda));

%--------------------------------------------------------------------------
% Reset random numbers generator
%--------------------------------------------------------------------------
seed = 1;
s = RandStream('mt19937ar','Seed', seed);
RandStream.setGlobalStream(s);


%% ------------------------------------------------------------------------
% Define Problem
%--------------------------------------------------------------------------

N_InnerOuter    = 500:500:3000;
N_minmax = 10000:10000:100000;


initialise_problem   = str2func('init_tc_so_unconstr_IAC_2018');
initialise_algorithm = str2func('init_algo_so_ebro_IAC_2018_interior_point_1pop'); 



% vector of objective functions: [TC1 TC2 TC3 ... TC14]



N_pop_algo = 1;

% constraint function TCC:
N_constraint = 1;


% Number of repeated evaluation (in order to evaluate the Success Rate)
N_iteration = 100;


% save results in a .txt file (0: don't; 1: do)
save_results = 1;



%% ------------------------------------------------------------------------
% RUN
%--------------------------------------------------------------------------
[problem] = initialise_problem();

[algo_minmax, algo_outer, algo_inner, algo_decomposition] = initialise_algorithm(problem);

for N_inner_i = N_InnerOuter
    
    
    algo_inner.par.nFeValMax = N_inner_i;
    algo_outer.par.nFeValMax = N_inner_i;
    
    
    for N_minmax_i = N_minmax
        
        algo_minmax.par_minmax.maxnfeval = N_minmax_i;
        
        %------------------------------------------------------------------
        % save results in a .txt file
        %------------------------------------------------------------------
        if save_results
            
            dir_save = '/RESULTS_minmax_IAC2018';
            mkdir(strcat(pwd, dir_save));
            
            name_txt = strcat('interior_point_1pop_f_minmax_IAC2018_nfevalInnerOuter',num2str(N_inner_i),'_nfevalTot',num2str(N_minmax_i));
            filename = fullfile(strcat(pwd, '/', dir_save,'/'), name_txt);
            [fileID_minmax, msg] = fopen(filename,'w');
            
            if fileID_minmax < 0
                error('Failed to open file "%s" because: "%s"', filename, msg);
            end
            
            str_archive = [];
            for len_arch = 1: sum([problem.dim_u   problem.dim_d  1])
                str_archive = [str_archive,' ', '%12.8f'];
            end
            str_archive = [str_archive, '\n'];
        end
        %------------------------------------------------------------------
        
        
        
        for k =1 : N_iteration
            
            
            [minmin, minmax, LIST, LIST_EXACT] = ebro(problem, algo_inner, algo_outer, algo_minmax, algo_decomposition);
            
            
            
            %==========================================================
            if save_results
                fprintf(fileID_minmax, str_archive, [minmax.f minmax.d minmax.u{1}]);
            end
            %==========================================================
            
        end
        
        %==============================================================
        if save_results
            fclose(fileID_minmax);
        end
        %==============================================================
    end
    
    
    
end



