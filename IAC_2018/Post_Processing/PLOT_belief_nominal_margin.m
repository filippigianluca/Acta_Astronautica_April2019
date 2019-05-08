% plot belief curves, minmax, nominal solutions
clear all; close all; clc

flag_belief_mass_1_sample      = 0;
flag_belief_mass_2_sample      = 1;
flag_belief_mass_deterministic = 0;
flag_belief_mass_arch_1        = 0;
flag_belief_mass_arch_2        = 1;
flag_belief_mass_arch_3        = 0;
flag_max_C_minmax              = 0;
flag_nominal_margin            = 1;



lineStyles=linspecer(length(nonzeros([flag_belief_mass_1_sample ...
    flag_belief_mass_2_sample ...
    flag_belief_mass_deterministic ...
    flag_belief_mass_arch_1 ...
    flag_belief_mass_arch_2 ...
    flag_belief_mass_arch_3 ...
    flag_max_C_minmax ])),1);


%% plot

ind_plot = 1;
ind_color = 1;

figure
hold on

load('u_nominal.mat');
load('d_minmax_nu400');

par.fix.time = 365;
par.fix.nu = 400;

% belief ENM for mass 1 sample
if flag_belief_mass_1_sample
    load('ebro_IAC_constr_SO_nfeval_decomposition_MASS1000_umaxMASS.mat','LIST')
    Belief_mass_1_sample = LIST;
    s(ind_plot) = stairs(Belief_mass_1_sample.F_Bel, Belief_mass_1_sample.Bel,'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','minmax 1 sample');
    
    ind_plot = ind_plot + 1;
    ind_color = ind_color + 1;
end

% belief ENM for mass 2 sample
if flag_belief_mass_2_sample
    load('ebro_IAC_constr_SO_nfeval_decomposition_2samples1000.mat','LIST')
    Belief_mass_2_sample = LIST;
    s(ind_plot) = stairs(Belief_mass_2_sample.F_Bel, Belief_mass_2_sample.Bel,'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','minmax 2 samples');
    
    ind_plot = ind_plot + 1;
    

    % margin solution for nominal deterministic approach
    if flag_nominal_margin
        f_nominal = IAC2018_only_mass2_for_belief(dminmax, u_nominal, par);
        s(ind_plot) = plot([f_nominal f_nominal], [0 1],'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','nominal d-minmax');

        ind_plot = ind_plot + 1;  
        
        output = IAC2018_obj_constr_u_nominal(dminmax, u_nominal, par);
        f_margin = output.f;
        s(ind_plot) = plot([f_margin f_margin], [0 1],'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','margin d-minmax');
        
        ind_plot = ind_plot + 1;
        ind_color = ind_color + 1;
    end
end



% belief ENM for mass 1 sample (optimal deterministic design)
if flag_belief_mass_deterministic
    load('ebro_IAC_constr_SO_decomposition_deterministic_design1000.mat','LIST','Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','opt det');
    s(ind_plot) = stairs(LIST.F_Bel, LIST.Bel);
    
    ind_plot = ind_plot + 1;
    

    % nominal-u
    f_nominal = IAC2018_only_mass2_for_belief(d_min_u_nominal, u_nominal, par);
    s(ind_plot) = plot([f_nominal f_nominal], [0 1],'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','nominal d-minmax');

    ind_plot = ind_plot + 1; 
    
    
    % margin solution for nominal deterministic approach
    if flag_nominal_margin
        par.fix.time = 365;
        par.fix.nu = 400;
        output = IAC2018_obj_constr_u_nominal(d_min_u_nominal, u_nominal, par);
        f_margin = output.f;
        s(ind_plot) = plot([f_margin f_margin], [0 1],'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','margin nominal');
        
        ind_plot = ind_plot + 1;
        ind_color = ind_color + 1;
    end
end


% belief ENM for mass 1 sample (1-archive Ad-Au)
if flag_belief_mass_arch_1
    load('ebro_IAC_constr_SO_decomposition_minmax_archive1.mat','LIST','Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','Arch 1');
    s(ind_plot) = stairs(LIST.F_Bel, LIST.Bel);
    
    ind_plot = ind_plot + 1;
    
    % nominal-u
    f_nominal = IAC2018_only_mass2_for_belief(d_minmax_non_convergence_1, u_nominal, par);
    s(ind_plot) = plot([f_nominal f_nominal], [0 1],'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','nominal d-minmax');

    ind_plot = ind_plot + 1; 
    
    
    % margin solution for 1-archive Ad-Au,
    if flag_nominal_margin
        load('d_u_minmax_non_convergence');
        par.fix.time = 365;
        par.fix.nu = 400;
        output = IAC2018_obj_constr_u_nominal(d_minmax_non_convergence_1, u_nominal, par);
        f_margin = output.f;
        s(ind_plot) = plot([f_margin f_margin], [0 1],'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','margin Arch 1');
        
        ind_plot = ind_plot + 1;
        ind_color = ind_color + 1;
    end
end


% belief ENM for mass 1 sample (2-archive Ad-Au)
if flag_belief_mass_arch_2
    load('ebro_IAC_constr_SO_decomposition_minmax_archive2.mat','LIST');
    load('d_u_minmax_non_convergence');
    s(ind_plot) = stairs(LIST.F_Bel, LIST.Bel,'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','Arch 2');
    
    ind_plot = ind_plot + 1;
    
    % nominal-u
    f_nominal = IAC2018_only_mass2_for_belief(d_minmax_non_convergence_2, u_nominal, par);
    s(ind_plot) = plot([f_nominal f_nominal], [0 1],'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','nominal d-minmax');

    ind_plot = ind_plot + 1; 
    
    
    % margin solution for 2-archive Ad-Au,
    if flag_nominal_margin
        
        par.fix.time = 365;
        par.fix.nu = 400;
        output = IAC2018_obj_constr_u_nominal(d_minmax_non_convergence_2, u_nominal, par);
        f_margin = output.f;
        s(ind_plot) = plot([f_margin f_margin], [0 1],'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','margin Arch 2');
        
        ind_plot = ind_plot + 1;
        ind_color = ind_color + 1;
    end
end


% belief ENM for mass 1 sample (3-archive Ad-Au)
if flag_belief_mass_arch_3
    load('ebro_IAC_constr_SO_decomposition_minmax_archive3.mat','LIST');
    load('d_u_minmax_non_convergence');
    s(ind_plot) = stairs(LIST.F_Bel, LIST.Bel,'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','Arch 3');
    
    ind_plot = ind_plot + 1;
    
    % nominal-u
    f_nominal = IAC2018_only_mass2_for_belief(d_minmax_non_convergence_3, u_nominal, par);
    s(ind_plot) = plot([f_nominal f_nominal], [0 1],'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','nominal d-minmax');

    ind_plot = ind_plot + 1; 
    
    
    % margin solution for 3-archive Ad-Au,
    if flag_nominal_margin
        
        par.fix.time = 365;
        par.fix.nu = 400;
        output = IAC2018_obj_constr_u_nominal(d_minmax_non_convergence_3, u_nominal, par);
        f_margin = output.f;
        s(ind_plot) = plot([f_margin f_margin], [0 1],'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','margin Arch 3');
        
        ind_plot = ind_plot + 1;
        ind_color = ind_color + 1;
    end
end


% u max violation constraint for d minmax
if flag_max_C_minmax
    load('results_functions_uminmax_umaxC.mat')
    s(ind_plot) = plot([mass_dminmax_umaxC mass_dminmax_umaxC], [0 1],'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','max C minmax');
    
    ind_plot = ind_plot + 1;
    ind_color = ind_color + 1;
end


%
if flag_max_C_minmax
    load('u_nominal.mat')
    load('minimum_d_nominal_u')
    s(ind_plot) = plot([f_min_u_nominal f_min_u_nominal], [0 1],'Color',lineStyles(ind_color,:),'linewidth',1,'DisplayName','max C');
    
    ind_plot = ind_plot + 1;
    ind_color = ind_color + 1;
end


%
% load('minimum_d_nominal_u.mat')
% plot([13.44 13.44], [0 1])
% ind_plot = ind_plot + 1;





lgd = legend(s(1,:));
hold off