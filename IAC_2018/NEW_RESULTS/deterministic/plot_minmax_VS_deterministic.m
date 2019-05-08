% nu= 600 Gbits

clear all; close all; clc

nu = 600;
par.fix.time = 365;
par.fix.nu = nu;


func_handle = @Acta_Astronautica_objconstr1;


% minmax solution
load('Acta_nu600_so_nfeval_minmax_3000000nfeval_inner_outer_20000', 'minmax');
minmax_d = minmax.d;
minmax_u = minmax.u{1};
minmax_f = minmax.f{1};


% optimal deterministic design d
load('u_maxc_nu600_dminmax','u_max_c')
load('min_deterministic','d_min_det');
load('u_nominal.mat','u_nominal')
u_nominal = [u_nominal(1:12) u_nominal(14:end)];


% constrained worst uncertain u for objective function mass
load('constr_worst_case_d_det');
u_constr_worst_mass = worst_case.u{1, 1};
f_constr_worst_mass = worst_case.f{1};


% UNconstrained worst uncertain u for objective function mass
load('unconstr_worst_case_d_det');
u_UNconstr_worst_mass = worst_case.u{1, 1};
f_UNconstr_worst_mass = worst_case.f{1};


% worst uncertain u for constraint function data volume
load('worst_case_c_ddet_')
u_maxC_ddet = worst_case.u{1, 1};
maxC = worst_case.f{1, 1};

lineStyles=linspecer(5,1);


%% plot
hold on

p = plot([600 600], [10.5 13.5],'k');

% constrained minmax solution
load('map_u_info');
u_max_c = map_affine(u_max_c, map_u_info);
output_minmax = func_handle(minmax_d, u_max_c, par);
minmax_c = nu - output_minmax.c;

s(1) = scatter(minmax_c, minmax_f, 'filled', 'MarkerFaceColor',lineStyles(2,:),'DisplayName',strcat('minmax'));
s(1).SizeData = 300;

% deterministic solution
output_det = func_handle(d_min_det, u_nominal, par);
det_f = output_det.f;
det_c = nu - output_det.c;
s(2) = scatter(det_c, det_f, 'filled', 'MarkerFaceColor',lineStyles(1,:),'DisplayName',strcat('opt det'));
s(2).SizeData = 300;
s(2).Marker = 's';

% % deterministic solution - constr max mass 
% constr_output_det = func_handle(d_min_det, u_constr_worst_mass, par);
% constr_det_f = constr_output_det.f;
% constr_det_c = nu - constr_output_det.c;
% s3 = scatter(constr_det_c, f_constr_worst_mass, 'filled', 'MarkerFaceColor',lineStyles(3,:),'DisplayName',strcat('d det, constr max mass'));


% deterministic solution - UNconstr max mass 
UNconstr_output_det = func_handle(d_min_det, u_UNconstr_worst_mass, par);
UNconstr_det_f = UNconstr_output_det.f;
UNconstr_det_c = nu - UNconstr_output_det.c;
s(3) = scatter(UNconstr_det_c, f_UNconstr_worst_mass, 'filled', 'MarkerFaceColor',lineStyles(3,:),'DisplayName',strcat('d det, max mass'));
s(3).SizeData = 300;
s(3).Marker = 'h';


% deterministic solution - max constraint violation
output_maxC = func_handle(d_min_det, u_maxC_ddet, par);
f_maxC = output_maxC.f;
s(4) = scatter(maxC, f_maxC, 'filled', 'MarkerFaceColor',lineStyles(4,:),'DisplayName',strcat('d det, max C'));
s(4).SizeData = 300;
s(4).Marker = 'p';




%%
grid on
xlabel('Data Volume')
ylabel('Mass')
lgd = legend(s(1,:));  