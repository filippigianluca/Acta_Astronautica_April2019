% run N random design and take the worst and best case
clear all; close all; clc

threshold = 600;
rip_tot = 1000;

% load deterministic case
load('min_deterministic','d_min_det');
load('worst_case_c_ddet_')
d_opt_det = d_min_det;
u_maxC_d_opt_det = worst_case.u{1, 1}  ;  


% load minmax case
load(strcat('Acta_nu',num2str(threshold),'_so_nfeval_minmax_3000000nfeval_inner_outer_20000'),'minmax','problem_minmax');
minmax_d = minmax.d;
minmax_u = minmax.u{1};


ind=1;
for k=[1 3 17 18 19]
load('dsave')
if k ==1
    d_min_det;
else
    d_opt_det = dsave(k,:);
end
% d_opt_det = problem_minmax.lb_d + rand(12,1).*(problem_minmax.ub_d-problem_minmax.lb_d);
% d_opt_det = d_opt_det';


%% 
[state_minmax,  s_minmax]   = histogram_minmax_1(rip_tot, minmax_d, minmax_u);

[b, bb] = histogram_nominal_1(rip_tot, d_opt_det, u_maxC_d_opt_det);
state_nominal{ind,:}=b;
s_nominal{ind,:}=bb;
ind=ind+1;
end


%%


figure
hold on

h1_nominal1 = hist(state_nominal{1},[0 1 2]);
h1_nominal1 = h1_nominal1/sum(h1_nominal1);

h1_nominal2 = hist(state_nominal{2},[0 1 2]);
h1_nominal2 = h1_nominal2/sum(h1_nominal2);

h1_nominal3 = hist(state_nominal{3},[0 1 2]);
h1_nominal3 = h1_nominal3/sum(h1_nominal3);

h1_nominal4 = hist(state_nominal{4},[0 1 2]);
h1_nominal4 = h1_nominal4/sum(h1_nominal4);

h1_nominal5 = hist(state_nominal{5},[0 1 2]);
h1_nominal5 = h1_nominal5/sum(h1_nominal5);

h1_minmax = hist(state_minmax,[0 1 2]);
h1_minmax = h1_minmax/sum(h1_minmax);

bar([0 1 2],[h1_nominal1;h1_nominal2;h1_nominal3;h1_nominal4;h1_nominal5;h1_minmax]')
%---------------------------------------
set(gca,'XTickLabel',{'0','','1','','2'});
legend('show');


figure
hold on

h2_nominal1 = hist(s_nominal{1},[0 1 2]);
h2_nominal1 = h2_nominal1/sum(h2_nominal1);

h2_nominal2 = hist(s_nominal{2},[0 1 2]);
h2_nominal2 = h2_nominal2/sum(h2_nominal2);

h2_nominal3 = hist(s_nominal{3},[0 1 2]);
h2_nominal3 = h2_nominal3/sum(h2_nominal3);

h2_nominal4 = hist(s_nominal{4},[0 1 2]);
h2_nominal4 = h2_nominal4/sum(h2_nominal4);

h2_nominal5 = hist(s_nominal{5},[0 1 2]);
h2_nominal5 = h2_nominal5/sum(h2_nominal5);

h2_minmax = hist(s_minmax,[0 1 2]);
h2_minmax = h2_minmax/sum(h2_minmax);

bar([0 1 2],[h2_nominal1;h2_nominal2;h2_nominal3;h2_nominal4;h2_nominal5;h2_minmax]')
set(gca,'XTickLabel',{'0','','1','','2'});
legend('show');


