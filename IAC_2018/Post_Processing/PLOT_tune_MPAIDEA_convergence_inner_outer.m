% This script studies the convergence of MP-AIDEA for IAC2018 test-case for
% the minmax problem. Inner and outer loops are considered separatdly: for
% different fixed values of design d the function is maximised; for
% different values of the uncertain vector u, minimisations are run.
% Different parameters for MP-AIDEA are considered:
% - number of agents in the population;
% - fmincon set: "sqp" or "interior-point";
% - maximum number of local restart;
clear all; close all; clc


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin user-configurable inputs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of fixed design (or uncertain vector) from "D_matrix" (or "U_matrix")
matrix_n_fix = 18;

% threshold for constrained function 
nu = 400; % fix

% first point to be plotted from the vector results
start_point_plot = 1;

% choos maximisation in the inner loop or minimisation in the outer loop
choose_loop = 'outer';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end user-configurable inputs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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

% results .txt in H drive
addpath(genpath('H:\My Documents\RESULTS_UTOPIAE_RELIABLES\ACTA_ASTRONAUTICA_IAC2018'));

% tool for plot in smart-o2c repository
addpath(genpath('C:\Users\yhb17181\Documents\MY_GITHUB\smart-o2c\Other_Tools'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end add folders to path %
%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% 

if strcmp(choose_loop,'inner')
    load('D_matrix');
    init_tc   = str2func(strcat('init_tc_so_constr_IAC_2018_tune_MPAIDEA_fixed_d'));
    IAC_mask_constr = @IAC_constr_tune_fix_d;
    fix_name = 'D';
elseif strcmp(choose_loop,'outer')
    load('U_matrix');
    init_tc   = str2func(strcat('init_tc_so_constr_IAC_2018_tune_MPAIDEA_fixed_u'));
    IAC_mask_constr = @IAC_constr_tune_fix_u;
    fix_name = 'U';
end

[ par ] = init_tc();
lineStyles=linspecer(4,1);








%--------------------------------------------------------------------------
n_agents = 10;

fmincon_set = 'sqp';
load(strcat('tune_constr1obj_',choose_loop,'_30000nfeval_',num2str(n_agents),'agents_',fix_name,num2str(matrix_n_fix),fmincon_set), 'output');
output_sqp_10agents = output;
Mass_sqp_10agents = output_sqp_10agents.memories_record;

fmincon_set = 'interior-point';
load(strcat('tune_constr1obj_',choose_loop,'_30000nfeval_',num2str(n_agents),'agents_',fix_name,num2str(matrix_n_fix),fmincon_set), 'output');
output_interior_10agents = output;
Mass_interior_10agents = output_interior_10agents.memories_record;
%--------------------------------------------------------------------------
n_agents = 20;

fmincon_set = 'sqp';
load(strcat('tune_constr1obj_',choose_loop,'_30000nfeval_',num2str(n_agents),'agents_',fix_name,num2str(matrix_n_fix),fmincon_set), 'output');
output_sqp_20agents = output;
Mass_sqp_20agents = output_sqp_20agents.memories_record;

fmincon_set = 'interior-point';
load(strcat('tune_constr1obj_',choose_loop,'_30000nfeval_',num2str(n_agents),'agents_',fix_name,num2str(matrix_n_fix),fmincon_set), 'output');
output_interior_20agents = output;
Mass_interior_20agents = output_interior_20agents.memories_record;
%--------------------------------------------------------------------------

if strcmp(choose_loop,'inner')
    d = D(matrix_n_fix,:);
    par.d = par.lb_d' + d.*par.ub_d';
elseif strcmp(choose_loop,'outer')
    u = U(matrix_n_fix,:);
    map_info = get_map_info(par.lb_u{1}, par.ub_u{1});
    par.u = map_affine(u, map_info);
end

for j=1:size(output_sqp_10agents.memories_record, 1)    
    
    u_sqp = output_sqp_10agents.memories_record(j,1:end-2);   
    CONSTR_sqp_10agents(j) = IAC_mask_constr(u_sqp, par);
end

for j=1:size(output_interior_10agents.memories_record, 1)    
 
    u_interior = output_interior_10agents.memories_record(j,1:end-2);
    CONSTR_interior_10agents(j) = IAC_mask_constr(u_interior, par);
end


for j=1:size(output_sqp_20agents.memories_record, 1)    

    u_sqp = output_sqp_20agents.memories_record(j,1:end-2);
    CONSTR_sqp_20agents(j) = IAC_mask_constr(u_sqp, par);
end

for j=1:size(output_interior_20agents.memories_record, 1)    

    u_interior = output_interior_20agents.memories_record(j,1:end-2);
    CONSTR_interior_20agents(j) = IAC_mask_constr(u_interior, par);
end

if strcmp(choose_loop,'inner')
    Mass_sqp_10agents = - Mass_sqp_10agents;
    Mass_interior_10agents = - Mass_interior_10agents;
    Mass_sqp_20agents = - Mass_sqp_20agents;
    Mass_interior_20agents = - Mass_interior_20agents;
end
%% plot 

ff = figure;


subplot(2,1,1);
hold on
p1 = plot(output_sqp_10agents.memories_record(start_point_plot:end,end),       Mass_sqp_10agents(start_point_plot:end,end-1),'Color',lineStyles(1,:));
p2 = plot(output_interior_10agents.memories_record(start_point_plot:end,end),  Mass_interior_10agents(start_point_plot:end,end-1),'Color',lineStyles(2,:));
p3 = plot(output_sqp_20agents.memories_record(start_point_plot:end,end),       Mass_sqp_20agents(start_point_plot:end,end-1),'Color',lineStyles(3,:));
p4 = plot(output_interior_20agents.memories_record(start_point_plot:end,end),  Mass_interior_20agents(start_point_plot:end,end-1),'Color',lineStyles(4,:));

s1 = scatter(output_sqp_10agents.memories_record(start_point_plot:end,end),       Mass_sqp_10agents(start_point_plot:end,end-1),'filled');
s2 = scatter(output_interior_10agents.memories_record(start_point_plot:end,end),  Mass_interior_10agents(start_point_plot:end,end-1),'filled');
s3 = scatter(output_sqp_20agents.memories_record(start_point_plot:end,end),       Mass_sqp_20agents(start_point_plot:end,end-1),'filled');
s4 = scatter(output_interior_20agents.memories_record(start_point_plot:end,end),  Mass_interior_20agents(start_point_plot:end,end-1),'filled');

s1.MarkerFaceColor = lineStyles(1,:);
s2.MarkerFaceColor = lineStyles(2,:);
s3.MarkerFaceColor = lineStyles(3,:);
s4.MarkerFaceColor = lineStyles(4,:);

xlabel('nfeval')
ylabel('Mass')
legend('sqp 10agents','interior 10agents', 'sqp 20agents','interior 20agents')
title(strcat('Mass', choose_loop,' fixed D',num2str(matrix_n_fix),'\nu = ',num2str(nu)))





subplot(2,1,2);
hold on
p5 = plot(output_sqp_10agents.memories_record(start_point_plot:end,end),      CONSTR_sqp_10agents(start_point_plot:end),'Color',lineStyles(1,:));
p6 = plot(output_interior_10agents.memories_record(start_point_plot:end,end), CONSTR_interior_10agents(start_point_plot:end),'Color',lineStyles(2,:));
p7 = plot(output_sqp_20agents.memories_record(start_point_plot:end,end),      CONSTR_sqp_20agents(start_point_plot:end),'Color',lineStyles(3,:));
p8 = plot(output_interior_20agents.memories_record(start_point_plot:end,end), CONSTR_interior_20agents(start_point_plot:end),'Color',lineStyles(4,:));

s5 = scatter(output_sqp_10agents.memories_record(start_point_plot:end,end),      CONSTR_sqp_10agents(start_point_plot:end),'filled');
s6 = scatter(output_interior_10agents.memories_record(start_point_plot:end,end), CONSTR_interior_10agents(start_point_plot:end),'filled');
s7 = scatter(output_sqp_20agents.memories_record(start_point_plot:end,end),      CONSTR_sqp_20agents(start_point_plot:end),'filled');
s8 = scatter(output_interior_20agents.memories_record(start_point_plot:end,end), CONSTR_interior_20agents(start_point_plot:end),'filled');

s5.MarkerFaceColor = lineStyles(1,:);
s6.MarkerFaceColor = lineStyles(2,:);
s7.MarkerFaceColor = lineStyles(3,:);
s8.MarkerFaceColor = lineStyles(4,:);

xlabel('nfeval')
ylabel('Constr')
legend('sqp 10agents','interior 10agents', 'sqp 20agents','interior 20agents')
title(strcat('Constr', choose_loop,' fixed D',num2str(matrix_n_fix),'\nu = ',num2str(nu)))

% saveas(ff,strcat('IAC2018_inner_convergence_D',num2str(matrix_n_fix),'.png'));


