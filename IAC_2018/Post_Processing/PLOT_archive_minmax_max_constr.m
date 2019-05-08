% plot comparison SO-MO minmax solutions

clear all; close all; clc


%% %%%%%%%%%%%%%%%
% set parameters %
%%%%%%%%%%%%%%%%%%

function_handle = @Acta_Astronautica_objconstr1; %@IAC2018_obj_constr;

plot_exact_maxima_result_minmax = 0;

nu = 400;

load_file = [];

% 1 -> 1 population,  LR CR and F fixed
% 2 -> 3 populations, LR CR and F adaptive
type_simulation = 1;


MO_flag = 0;

unconstr_SO_flag = 1;

objfun_flag = 1;

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


figure


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
%     ind = find(nu==j);
    
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
    
    %     [fminmax, pos_minmax] = min(f_minmax(max_C_violation));
    f_minmax_feasible     = f_minmax(max_C_violation);
    d_minmax_feasible     = d_minmax(max_C_violation);
    u_minmax_feasible     = u_minmax(max_C_violation);
    u_max_constr_feasible = u_max_constr(max_C_violation);
    
    lineStyles=linspecer(length(max_C_violation)+1,1); 
    
    ind_colour = 1;
    ind_plot   = 1;
    aaa = find(max_C_violation==1);
    for iii = 1 : sum(max_C_violation)
        
        if objfun_flag == 1
            outpu_u   = function_handle(d_minmax_feasible{(iii)},  u_minmax_feasible{(iii)}, par);
            MASS = outpu_u.f;
            RES = - outpu_u.c + j;
            outpu_uC   = function_handle(d_minmax_feasible{(iii)}, u_max_constr_feasible{(iii)}, par);
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
        
        
        
        s(ind_plot) = scatter(RES, MASS,'filled','DisplayName',strcat('constr SO \nu = ',num2str(j)));
        s(ind_plot).MarkerFaceColor = lineStyles(ind_colour,:);
        s(ind_plot).SizeData = 100;
        
        s(ind_plot+1) = scatter(Max_RES, MASS_max_RES,'filled','DisplayName',strcat('constr SO max violation'));
        s(ind_plot+1).MarkerFaceColor = lineStyles(end,:);
        s(ind_plot+1).SizeData = 100;
        %    plot([j j], [0 80],'k'); %  plot([j j], [9.5 MASS],'k');
        plot([RES Max_RES], [MASS MASS_max_RES],'Color',lineStyles(end,:),'linewidth',1);
        
        
        ind_colour = ind_colour + 1;
        ind_plot   = ind_plot + 2;
    end
end


if unconstr_SO_flag
    load('IAC2018_3pop_1e6_2e4_2e4_sqp_feval1100000','ARCHIVE');
    par.fix.nu = 400;
    SO_result = function_handle(ARCHIVE{1, 1}{11, 1}  , ARCHIVE{1, 1}{11, 2}  , par  );
    s(end+1) = scatter(400-SO_result.c, SO_result.f,'DisplayName',strcat('unconstr SO' ));
    s(end).MarkerFaceColor = [0 0 0];
    s(end).MarkerEdgeColor = [0 0 0];
    s(end).SizeData = 100;
    s(end).Marker = 'p';
end


xlabel('E(V)')
ylabel('Mass')
lgd = legend(s(1,:));
title(strcat('constrained Mass minmax'))
