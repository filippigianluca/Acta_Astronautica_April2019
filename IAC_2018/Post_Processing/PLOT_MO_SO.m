% plot comparison SO-MO minmax solutions

clear all; close all; clc


%% %%%%%%%%%%%%%%%
% set parameters %
%%%%%%%%%%%%%%%%%%

plot_exact_maxima_result_minmax = 0;

nu = [380 390 400 410];

load_file = [];

% 1 -> 1 population,  LR CR and F fixed
% 2 -> 3 populations, LR CR and F adaptive
type_simulation = 1;


MO_flag          = 0;
unconstr_SO_flag = 0;
objfun_flag      = 1;
plot_flag        = 2; % 2: plot(min_E(V), max(Mass));
plot_number_flag = 1:4;%3;
nominal_flag     = 0;
%% %%%%%%%%%%%%%%%%%%%
% end set parameters %
%%%%%%%%%%%%%%%%%%%%%%





%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% start add folders to path %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all')
%warning('on','all')


% folders in UTOPIAR_RELIABLES repository
CURRENTPATH=pwd;
if isunix
    idcs   = strfind(CURRENTPATH,'/');
else
    idcs   = strfind(CURRENTPATH,'\');
end
ALLdir    = CURRENTPATH(1:idcs(end-3)-1);
MY_GITHUB = CURRENTPATH(1:idcs(end-2)-1);
if isunix
    Stuff_Danda = strcat(MY_GITHUB,'/utopiae-reliables/Stuff_DANDA');
    IAC2018     = strcat(MY_GITHUB,'/utopiae-reliables/IAC_2018');
else
    Stuff_Danda = strcat(MY_GITHUB,'\utopiae-reliables\Stuff_DANDA');
    IAC2018     = strcat(MY_GITHUB,'\utopiae-reliables\IAC_2018');
end

rmpath(genpath(ALLdir));
addpath(genpath(IAC2018));
addpath(genpath(Stuff_Danda));


% % results .txt in H drive
% addpath(genpath('H:\My Documents\RESULTS_UTOPIAE_RELIABLES\ACTA_ASTRONAUTICA_IAC2018'));
%
% % tool for plot in smart-o2c repository
addpath(genpath('C:\Users\yhb17181\Documents\MY_GITHUB\smart-o2c\Other_Tools'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end add folders to path %
%%%%%%%%%%%%%%%%%%%%%%%%%%%


if plot_flag ==1
    lineStyles=linspecer(length(nu)+1,1);
elseif plot_flag == 2
    lineStyles=linspecer(length(nu),1);
end


%% MO
MO_p = 0;
if MO_flag
    load('IAC_2018_mo_nfeval_minmax_1000000nfeval_inner_outer_20000', 'fminmax')
    s(1) = scatter(-fminmax(:,2), fminmax(:,1),'filled','DisplayName',strcat('MO'));
    s(1).MarkerFaceColor = lineStyles(1,:);
    MO_p = 1;
end

%% SO
for j = nu
    ind = find(nu==j);
    
    hold on
    
    % load('nu390_20000inout_cnstr_1obj_nfeval820015','ARCHIVE');
    
    if isempty(load_file)
        %     if j==390
        %         load(strcat('nu',num2str(j),'_IAC_2018_so_constr_nfeval_minmax_2000000nfeval_inner_outer_20000'),'ARCHIVE');
        %         Arch = ARCHIVE{1};
        %     else
        %         load(strcat('nu',num2str(j),'_IAC_2018_so_constr_nfeval_minmax_2000000nfeval_inner_outer_20000'),'minmax');
        %         Arch = minmax.output.ARCHIVE{1, 1};
        %     end
        if type_simulation==1
            load(strcat('nu',num2str(j),'_20000inout_cnstr_1pop_1obj_LR_F_CR_fixed'),'ARCHIVE');
            Arch = ARCHIVE{1};
        elseif type_simulation==2
            load(strcat('nu',num2str(j),'_20000inout_cnstr_3pop_1obj_LR_F_CR_adaptive'),'ARCHIVE');
            Arch = ARCHIVE{1, 1};
        end
    else
        load(load_file);
        Arch = ARCHIVE{1};
    end
    
    
    
    for i=2:size(Arch, 1)
        par.fix.nu = j;
        par.fix.time = 365;
        f_minmax(i-1) = Arch{i,4};
        max_C_violation(i-1) = Arch{i,6}<=0;
        N_minmax(i-1) = Arch{i,7};
        
        d_minmax{i-1} = Arch{i,1};
        u_minmax{i-1} = Arch{i,2};
        u_max_constr{i-1} = Arch{i,3};
        
    end
    
    [fminmax, pos_minmax] = min(f_minmax(max_C_violation));
    f_minmax_feasible     = f_minmax(max_C_violation);
    d_minmax_feasible     = d_minmax(max_C_violation);
    u_minmax_feasible     = u_minmax(max_C_violation);
    u_max_constr_feasible = u_max_constr(max_C_violation);
    
    if objfun_flag == 1
        outpu_u   = IAC2018_obj_constr(d_minmax_feasible{pos_minmax},  u_minmax_feasible{pos_minmax}, par);
        MASS = outpu_u.f;
        RES = - outpu_u.c + j;
        outpu_uC   = IAC2018_obj_constr(d_minmax_feasible{pos_minmax}, u_max_constr_feasible{pos_minmax}, par);
        MASS_max_RES = outpu_uC.f;
        Max_RES = - outpu_uC.c + j;
    else
        MASS   = CUBESAT_5subsystems_MASS(Arch{1+pos_minmax,1}, Arch{1+pos_minmax,2}, par);
        Constr = constr_mask_CUBESAT_5subsystems_RES(Arch{1+pos_minmax,1}, Arch{1+pos_minmax,2}, par);
        %     Max_Constr = constr_mask_CUBESAT_5subsystems_RES(d_minmax, u_max_constr, par);
        MASS_max_RES   = CUBESAT_5subsystems_MASS(Arch{1+pos_minmax,1}, Arch{1+pos_minmax,3}, par);
        RES     = - CUBESAT_5subsystems_RES(Arch{1+pos_minmax,1}, Arch{1+pos_minmax,2}, par);
        Max_RES = - CUBESAT_5subsystems_RES(Arch{1+pos_minmax,1}, Arch{1+pos_minmax,3}, par);
    end
    
    if any(ind == plot_number_flag)
        plot([j j], [8 14],'k'); % plot([j j], [10 80],'k'); %  plot([j j], [9.5 MASS],'k');
        if plot_flag == 1
            
            plot([RES Max_RES], [MASS MASS_max_RES],'Color',lineStyles(end,:),'linewidth',1);
            s(ind+MO_p) = scatter(RES, MASS,'filled','DisplayName',strcat('constr SO \nu = ',num2str(j)));
            s(ind+MO_p).MarkerFaceColor = lineStyles(ind,:);
            s(ind+MO_p).SizeData = 200;
            
            s(ind+MO_p+1) = scatter(Max_RES, MASS_max_RES,'filled','DisplayName',strcat('constr SO max violation'));
            s(ind+MO_p+1).MarkerFaceColor = lineStyles(ind,:);
            s(ind+MO_p+1).SizeData = 200;
            s(ind+MO_p+1).Marker = 'd';
            
           
            
        elseif plot_flag == 2
            s(ind+MO_p) = scatter(Max_RES, MASS,'filled','DisplayName',strcat('constr SO \nu = ',num2str(j)));
            s(ind+MO_p).MarkerFaceColor = lineStyles(ind,:);
            s(ind+MO_p).SizeData = 200;
        end
        
    end
% dminmax=d_minmax_feasible{pos_minmax};
% uminmax=u_minmax_feasible{pos_minmax};
% umaxC_dminmax=u_max_constr_feasible{pos_minmax};
% save(strcat('minmax_nu',num2str(j)),'dminmax','uminmax','umaxC_dminmax')
end


if unconstr_SO_flag
    load('IAC2018_3pop_1e6_2e4_2e4_sqp_feval1100000','ARCHIVE');
    par.fix.nu = 400;
    d_unconstr_minmax = ARCHIVE{1, 1}{11, 1};
    u_unconstr_minmax = ARCHIVE{1, 1}{11, 2};
    SO_result = IAC2018_obj_constr(d_unconstr_minmax, u_unconstr_minmax, par  );
    s(end+1) = scatter(400-SO_result.c, SO_result.f,'DisplayName',strcat('unconstr SO' ));
    s(end).MarkerFaceColor = [0 0 0];
    s(end).MarkerEdgeColor = [0 0 0];
    s(end).SizeData = 100;
    s(end).Marker = 'p';
       
    % worst case constraint for unconstrained-d solution
    load('max_C_d_unconstrained_40e4.mat')
    d_unconstr_max_C = max.d;
    u_unconstr_max_C = max.u{1};
    SO_result_max_C = IAC2018_obj_constr(d_unconstr_max_C, u_unconstr_max_C, par  );
    s(end+1) = scatter(400-SO_result_max_C.c, SO_result_max_C.f,'DisplayName',strcat('min E(V) unconstr SO' ));
    
    % worst case mass for unconstrained-d solution
    load('max_Mass_d_unconstrained_40e4.mat')
    SO_result_max_M = IAC2018_obj_constr(max.d, max.u{1}, par  );
    s(end+1) = scatter(400-SO_result_max_M.c, SO_result_max_M.f,'DisplayName',strcat('min E(V) unconstr SO' ));
end


if nominal_flag
    load('u_nominal')
    u_nominal;                                                   % nominal u
    load('minimum_d_nominal_u')
    d_opt_unominal    = d_min_u_nominal;                         % optimal d with fixed nominal u
    Mass_opt_unominal = f_min_u_nominal;                         % Mass - nominal u, optimal d
    output = IAC2018_obj_constr(d_opt_unominal, u_nominal, par);
    RES_mass_opt_unominal = 400 - output.c;                      % DV   - nominal u, optimal d
    load('ebro_IAC_max_CONSTR_dmin_nominal_u')
    max_constr = worst_case;
    u_max_RES = worst_case.u{1};                                 % u for min DV with optimal d with nominal u
    RES_max_d_opt_unominal = -worst_case.f{1};                   % DV    - u for min DV with optimal d with nominal u
    output = IAC2018_obj_constr(d_opt_unominal, u_max_RES, par);
    Mass_res_max_d_opt_unominal = output.f;                      % Mass  - u for min DV with optimal d with nominal u
    load('ebro_IAC_max_MASS_dmin_nominal_u')
    u_max_Mass = worst_case.u{1, 1};                             % u for max Mass with optimal d with nominal u
    output = IAC2018_obj_constr(d_opt_unominal, u_max_Mass, par);
    Mass_max_d_opt_unominal = output.f;                          % Mass  - u for max Mass with optimal d with nominal u
    RES_Massmax_d_opt_unominal = 400 -output.c;                  % DV    - u for max Mass with optimal d with nominal u
    
    s(end+1) = scatter(RES_mass_opt_unominal,  Mass_opt_unominal,'DisplayName',strcat('optimal design with nominal u' ));
    s(end+1) = scatter(RES_max_d_opt_unominal, Mass_res_max_d_opt_unominal,'DisplayName',strcat('worst DV with optimal d with nominal u' ));
    s(end+1) = scatter(RES_Massmax_d_opt_unominal, Mass_max_d_opt_unominal,'DisplayName',strcat('worst Mass with optimal d with nominal u' ));
end


xlabel('E(V)')
ylabel('Mass')
lgd = legend(s(1,:));
title(strcat('constrained Mass minmax'))
