% plot Delta F in for each parameter
clear all; close all; clc

rmpath(genpath('C:\Users\yhb17181\Documents\LOCAL_REPOSITORY_SMARTO2C\smart-o2c\Problems\EBRO\WorkInProgress'))
addpath(genpath('C:\Users\yhb17181\Documents\LOCAL_REPOSITORY_SMARTO2C\smart-o2c\Problems\EBRO\WorkInProgress\EXACT_CUBESAT'))

[lb_d, ub_d, lb_u, ub_u, par] = bound_SPACECRAFT();

d_input = [];

N_d = 50;   % number of design samples
N_u = 100;   % number of samples in u_i interval



system = 'all';
objective_fun = 'mass';
% plot_figure = 0 --> no plot
% plot_figure = 1 --> plot
plot_figure     = 1;
plot_figure_std = 0;

if strcmp(system, 'aocs')
%--------------------------------------------------------------------------
% AOCS
func = @space_aocs;
[lb_d, ub_d, lb_u, ub_u] = bound_AOCS();
%--------------------------------------------------------------------------
elseif strcmp(system, 'ttc')
% TTC
func = @space_ttc;
[lb_d, ub_d, lb_u, ub_u] = bound_TTC();
%--------------------------------------------------------------------------
elseif strcmp(system, 'power')
% POWER
func = @space_power;
[lb_d, ub_d, lb_u, ub_u] = bound_POWER();
%--------------------------------------------------------------------------
elseif strcmp(system, 'all')
% ALL SYSTEM
func = @mask_spacecraft;
[lb_d, ub_d, lb_u, ub_u, par] = bound_SPACECRAFT();
lb = [lb_d lb_u];
ub = [ub_d ub_u];
par.dim_d = length(lb_d);
%--------------------------------------------------------------------------
end



d_input = ub_d;
d_input([1 6 9 12 13 14 16 18 22 24 26 28 29 30 31]) = [40,0.000500000000000000,200,30.1848785240731,11.2177542853998,-59.0112904291059,0.384940165924359,5,1,0.861044018717549,0.105936414254059,0.946947460004865,3.04702345973368,0.400000000000000,2.50000000000000];


for j=[1 5 6]%:length(lb_u)

    if plot_figure == 1
        figure
    end
    
    if ~isempty(d_input)
        N_d = 1;
    end
    
    for i=1:N_d
        
        if ~isempty(d_input)
            d = d_input;
        else
            d = lb_d + rand(1,length(lb_d)).*(ub_d-lb_d); 
        end
        
        
        for ind_ran = 1:20
            
        u = lb_u + rand(1,length(lb_u)).*(ub_u-lb_u);
        
        for k=1:N_u+1
            u(j) = lb_u(j) + (k-1)/N_u*(ub_u(j)-lb_u(j));
            if strcmp(system, 'all')
                x = [d, u];
                [M] = func(x, par);
            else
                [M,P,info] = func(d, u);
            end
            
            if strcmp (objective_fun, 'mass') == 1
                F(i, j, k) = M; 
                F_plot{ind_ran}(i,k) = M;
            elseif strcmp (objective_fun, 'power') == 1
                F(i, j, k) = P; 
                F_plot(i,k) = P;
            end
              
        end
        
        end
        

        %----------------------------
        standard_deviation(j,i) = std(F(i, j,:));
        monotone(j,i) = issorted(F(i, j,:)) || issorted(-F(i, j,:));
        minimum(j,i) = min(F(i, j,:));
        maximum(j,i) = max(F(i, j,:));
%         figure
%         errorbar(1:length(minimum),(minimum+maximum)/2,maximum - minimum)
        %----------------------------
        
        if plot_figure == 1
            for ind_ran = 1:10
            U(i,:) = u;
            D(i,:) = d;
            hold on
            plot(F_plot{ind_ran}(i,:));
            end
        end
    end
    standard_deviation_global(j,1) = mean(standard_deviation(j,:));
    monotone_global(j,1) = mean(monotone(j,:));   

    
        if plot_figure_std == 1
            clear U D
            title(['std:',num2str(standard_deviation_global(j,1))])
            drawnow
        end
end