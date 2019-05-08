% plot for different thresholds in the constraint the minmax solution: 
% 1. d-minmax, u-max-constraint


% clear all; close all; clc

%% defin parameters

nu_vec   = 500;
N_in_out = 10000;
N_tot    = 2000000;



%% run
for nu = nu_vec
    
    clearvars -except nu N_tot N_in_out
    
    
    load(strcat('Acta_nu',num2str(nu),'_so_nfeval_minmax_',num2str(N_tot),'nfeval_inner_outer_',num2str(N_in_out)))
    
    flag_plot_mpaidea = 0;
    
    lineStyles=linspecer(5,1);
    
    f_minmax_convergence = cell2mat(minmax.archiveG(2:end,4));
    f_minmax_convergence_cc = cell2mat(minmax.archiveG_cc(2:end,4));
    c_max_convergence = nu-cell2mat(minmax.archiveG_cc(2:end,6));
    N = cell2mat(minmax.archiveG(2:end,7));
    
    
    % first point
    f_uf_convergence1 = cell2mat(minmax.archiveG(2,9));
    N_uf_convergence1 = 1 : N(1)/length(f_uf_convergence1): N(1);
    
    hold on
    % s1_x =[N(1)-length(f_uf_convergence1):N(1)-1];
    
    yyaxis left
    plot(N(2:end),f_minmax_convergence_cc(2:end))
    if flag_plot_mpaidea
        p11 = plot(N_uf_convergence1, f_uf_convergence1);%,'Color',lineStyles(1,:),'linewidth',2);
    end
    % s   = scatter(N(1), f_minmax_convergence(1), 'filled');%, 'MarkerFaceColor',lineStyles(4,:), 'SizeData', 80);
    % s_cc = scatter(N(1), f_minmax_convergence_cc(1), 'filled');%, 'MarkerFaceColor',lineStyles(3,:));
    
    
    yyaxis right
    plot(N(2:end),c_max_convergence(2:end))
    % s_c = scatter(N(1), c_max_convergence(1), 'filled');%, 'MarkerFaceColor',lineStyles(5,:));
    
    N_cv =1;
    N_iter  = 1;
    for i=3:size(minmax.archiveG,1)
        
        f_df_convergence = minmax.archiveG{i,8};
        f_uf_convergence = minmax.archiveG{i,9};
        
        N_df_convergence =  N(i-2)  +1              : ((N(i-2)+N(i-1))/2-N(i-2))/length(f_df_convergence)  : N(i-2)+(N(i-1)-N(i-2))/2;
        N_uf_convergence = N(i-2)+(N(i-1)-N(i-2))/2 +1: ((N(i-1)-N(i-2))/2)/length(f_uf_convergence)         : N(i-1) ;
        N_df_convergence(end) = N_uf_convergence(1);
        N_uf_convergence(end) = N(i-1);
        
        yyaxis left
        if flag_plot_mpaidea
            p22(N_iter) = plot(N_df_convergence, f_df_convergence);%,'Color',lineStyles(2,:),'linewidth',2);
            p11(N_iter) = plot(N_uf_convergence, f_uf_convergence);%,'Color',lineStyles(1,:),'linewidth',2);
        end
        
        s_f(N_iter) = scatter(N(i-1), f_minmax_convergence_cc(i-1), 'filled');%, 'MarkerFaceColor',lineStyles(4,:), 'SizeData', 80);
        s_cc(N_iter) = scatter(N(i-1), f_minmax_convergence_cc(i-1), 'filled');%, 'MarkerFaceColor',lineStyles(3,:));
        
        yyaxis right
        s_c(N_iter) = scatter(N(i-1), c_max_convergence(i-1), 'filled');%, 'MarkerFaceColor',lineStyles(3,:));
        if c_max_convergence(i-1)<nu
            s_cv(N_iter) = scatter(N(i-1), c_max_convergence(i-1), 'filled','r');
            s_cv(N_iter).SizeData=100; %(N_cv)
            s_cv(N_iter).Marker = 's';
        end
            N_iter=N_iter+1;
    end
    
    
end



grid on

yyaxis left
ylabel('mass')

yyaxis right
ylabel('data volume')