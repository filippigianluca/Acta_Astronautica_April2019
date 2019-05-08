% plot belief curves, minmax, nominal solutions
clear all; close all; clc


function_handle_margin = @margin_Acta_Astronautica_objconstr1;
nu = 600;
Nsamples_minmax = 3;
Nsamples_Arch = 2;

load('min_deterministic')
darchive = d_min_det;
load('u_nominal_20')


%% plot

lineStyles=linspecer(2,1);


ind_plot = 1;
ind_color = 1;

figure
hold on

load('u_nominal.mat');
load('d_minmax_nu400');

par.fix.time = 365;
par.fix.nu = nu;


% belief curve minmax
load(strcat(num2str(Nsamples_minmax),'sample1000'),'LIST')
Belief_mass_minmax = LIST;
s(1) = stairs(Belief_mass_minmax.F_Bel, Belief_mass_minmax.Bel,'Color',lineStyles(1,:),'linewidth',2,'DisplayName','minmax');


% belief curve archive
load(strcat('Archive_',num2str(Nsamples_Arch),'sample1000'),'LIST')
Belief_mass_Arch = LIST;
s(2) = stairs(Belief_mass_Arch.F_Bel, Belief_mass_Arch.Bel,'Color',lineStyles(2,:),'linewidth',2,'DisplayName','archive');


% margin solution for nominal deterministic approach
f_margin = function_handle_margin(darchive, u_nominal_20, par);
s(3) = plot([f_margin(2) f_margin(2)], [0 1],'Color',lineStyles(2,:),'linewidth',2,'LineStyle','-.','DisplayName','nominal d-minmax');
s(4) = plot([f_margin(1) f_margin(1)], [0 1],'Color',lineStyles(2,:),'linewidth',2,'DisplayName','margin d-minmax');


grid on
lgd = legend(s(1,:));
hold off