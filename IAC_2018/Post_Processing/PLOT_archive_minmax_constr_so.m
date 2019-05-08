% plot results of minmax problem:
clear all; clc


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin user-configurable inputs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of fixed vector from "D_matrix"
D_ni = 9;

% threshold for constrained function 
nu = 390;

% first point to be plotted from the vector results
start_point_plot = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end user-configurable inputs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
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






%% run

init_tc   = str2func(strcat('init_tc_so_constr_IAC_2018'));
[ par ] = init_tc();
lineStyles=linspecer(4,1);


%--------------------------------------------------------------------------
% last results: 1-/02/2019
load('nu400_20000inout_cnstr_1obj_nfeval3020088','ARCHIVE')
% load('nu400_1000inout_cnstr_1obj_nfeval1493123','ARCHIVE')
% load('nu400_20000inout_cnstr_1obj_nfeval1420021','ARCHIVE')
% load('nu400__10000inout_cnstr_1obj_nfeval1480041','ARCHIVE')
% load('nu400__5000inout_cnstr_1obj_nfeval1475075','ARCHIVE')
% load('nu410_20000inout_cnstr_1obj_nfeval1420021','ARCHIVE')

Arch = ARCHIVE{1, 1};

%--------------------------------------------------------------------------
par.fix.nu = nu;
for i=2:size(Arch, 1)
    f_minmax(i-1) = Arch{i,4};
    N_minmax(i-1) = Arch{i,7};
    
    d_minmax = Arch{i,1};
    u_minmax = Arch{i,2};
    u_max_constr = Arch{i,3};
    
    MASS(i-1) = CUBESAT_5subsystems_MASS(d_minmax, u_minmax, par);
    Constr(i-1) = constr_mask_CUBESAT_5subsystems_RES(d_minmax, u_minmax, par);
    Max_Constr(i-1) = constr_mask_CUBESAT_5subsystems_RES(d_minmax, u_max_constr, par);
    RES(i-1) = -CUBESAT_5subsystems_RES(d_minmax, u_minmax, par);
    Max_RES(i-1) = -CUBESAT_5subsystems_RES(d_minmax, u_max_constr, par);
end



%% plot 

figure


subplot(2,1,1);
hold on
p1 = plot(output_sqp_10agents.memories_record(start_point_plot:end,end), -output_sqp_10agents.memories_record(start_point_plot:end,end-1),'Color',lineStyles(1,:));
p2 = plot(output_interior_10agents.memories_record(start_point_plot:end,end), -output_interior_10agents.memories_record(start_point_plot:end,end-1),'Color',lineStyles(2,:));
p3 = plot(output_sqp_20agents.memories_record(start_point_plot:end,end), -output_sqp_20agents.memories_record(start_point_plot:end,end-1),'Color',lineStyles(3,:));
p4 = plot(output_interior_20agents.memories_record(start_point_plot:end,end), -output_interior_20agents.memories_record(start_point_plot:end,end-1),'Color',lineStyles(4,:));

s1 = scatter(output_sqp_10agents.memories_record(start_point_plot:end,end), -output_sqp_10agents.memories_record(start_point_plot:end,end-1),'filled');
s2 = scatter(output_interior_10agents.memories_record(start_point_plot:end,end), -output_interior_10agents.memories_record(start_point_plot:end,end-1),'filled');
s3 = scatter(output_sqp_20agents.memories_record(start_point_plot:end,end), -output_sqp_20agents.memories_record(start_point_plot:end,end-1),'filled');
s4 = scatter(output_interior_20agents.memories_record(start_point_plot:end,end), -output_interior_20agents.memories_record(start_point_plot:end,end-1),'filled');

s1.MarkerFaceColor = lineStyles(1,:);
s2.MarkerFaceColor = lineStyles(2,:);
s3.MarkerFaceColor = lineStyles(3,:);
s4.MarkerFaceColor = lineStyles(4,:);

xlabel('nfeval')
ylabel('Mass')
legend('sqp 10agents','interior 10agents', 'sqp 20agents','interior 20agents')
title(strcat('Mass, INNER fixed D',num2str(D_ni),'\nu = ',num2str(nu)))





subplot(2,1,2);
hold on
p5 = plot(output_sqp_10agents.memories_record(start_point_plot:end,end), CONSTR_sqp_10agents(start_point_plot:end),'Color',lineStyles(1,:));
p6 = plot(output_interior_10agents.memories_record(start_point_plot:end,end), CONSTR_interior_10agents(start_point_plot:end),'Color',lineStyles(2,:));
p7 = plot(output_sqp_20agents.memories_record(start_point_plot:end,end), CONSTR_sqp_20agents(start_point_plot:end),'Color',lineStyles(3,:));
p8 = plot(output_interior_20agents.memories_record(start_point_plot:end,end), CONSTR_interior_20agents(start_point_plot:end),'Color',lineStyles(4,:));

s5 = scatter(output_sqp_10agents.memories_record(start_point_plot:end,end), CONSTR_sqp_10agents(start_point_plot:end),'filled');
s6 = scatter(output_interior_10agents.memories_record(start_point_plot:end,end), CONSTR_interior_10agents(start_point_plot:end),'filled');
s7 = scatter(output_sqp_20agents.memories_record(start_point_plot:end,end), CONSTR_sqp_20agents(start_point_plot:end),'filled');
s8 = scatter(output_interior_20agents.memories_record(start_point_plot:end,end), CONSTR_interior_20agents(start_point_plot:end),'filled');

s5.MarkerFaceColor = lineStyles(1,:);
s6.MarkerFaceColor = lineStyles(2,:);
s7.MarkerFaceColor = lineStyles(3,:);
s8.MarkerFaceColor = lineStyles(4,:);

xlabel('nfeval')
ylabel('Constr')
legend('sqp 10agents','interior 10agents', 'sqp 20agents','interior 20agents')
title(strcat('Constr, INNER fixed D',num2str(D_ni),'\nu = ',num2str(nu)))

