% run N random design and take the worst and best case


threshold = 600;
rip_tot = 5000;

% load deterministic case
load('min_deterministic','d_min_det');
load('worst_case_c_ddet_')
d_opt_det = d_min_det;
u_maxC_d_opt_det = worst_case.u{1, 1}  ;  


% load minmax case
load(strcat('Acta_nu',num2str(threshold),'_so_nfeval_minmax_3000000nfeval_inner_outer_20000'),'minmax','problem_minmax');
minmax_d = minmax.d;
minmax_u = minmax.u{1};


state_0_min = inf; state_1_min = inf; state_2_min = inf;
state_0_max = 0;   state_1_max = 0;   state_2_max = 0;
s_0_min     = inf; s_1_min     = inf; s_2_min     = inf;
s_0_max     = 0;   s_1_max     = 0;   s_2_max     = 0;


for k=1:300

d_opt_det = problem_minmax.lb_d + rand(12,1).*(problem_minmax.ub_d-problem_minmax.lb_d);
d_opt_det = d_opt_det';

% dsave(k,:)=d_opt_det;
%% 
[state_minmax,  s_minmax]   = histogram_minmax_1(rip_tot, minmax_d, minmax_u);

[state_nominal, s_nominal] = histogram_nominal_1(rip_tot, d_opt_det, u_maxC_d_opt_det);

a=length(state_nominal(state_nominal==0));
b=length(state_nominal(state_nominal==1));
c=length(state_nominal(state_nominal==2));

aa=length(s_nominal(s_nominal==0));
bb=length(s_nominal(s_nominal==1));
cc=length(s_nominal(s_nominal==2));

if a < state_0_min; state_0_min = a; end
if b < state_1_min; state_1_min = b; end
if c < state_2_min; state_2_min = c; end

if a > state_0_max; state_0_max = aa; end
if b > state_1_max; state_1_max = bb; end
if c > state_2_max; state_2_max = cc; end



if aa < s_0_min; s_0_min = a; end
if bb < s_1_min; s_1_min = b; end
if cc < s_2_min; s_2_min = c; end

if aa > s_0_max; s_0_max = aa; end
if bb > s_1_max; s_1_max = bb; end
if cc > s_2_max; s_2_max = cc; end
end




%%
figure
hold on

h1_nominal_1 = hist([zeros(1,state_0_max) ones(1,state_1_max) 2*ones(1,state_2_max)],[0 1 2]);
h1_nominal_1 = h1_nominal_1/sum(h1_nominal_1);

h1_nominal_2 = hist([zeros(1,state_0_min) ones(1,state_1_min) 2*ones(1,state_2_min)],[0 1 2]);
h1_nominal_2 = h1_nominal_2/sum(h1_nominal_2);

h1_nominal = [h1_nominal_1; h1_nominal_2];

h1_minmax = hist(state_minmax,[0 1 2]);
h1_minmax = h1_minmax/sum(h1_minmax);

% bar([0 1 2],[max(h1_nominal);h1_minmax]')
% bar([0 1 2],[min(h1_nominal);h1_minmax]')
bar([0 1 2],[h1_nominal_1;h1_nominal_2;h1_minmax]')
% bar([0 1 2],[h1_nominal_2;h1_minmax]')
%---------------------------------------

legend('show');




%%
figure
hold on

h2_minmax = hist(s_minmax,[0 1 2]);
h2_minmax = h2_minmax/sum(h2_minmax);
% bar(h1, 'DisplayName', 'N transitions');

% legend('show');
% title('min-max')



%%
h2_nominal_1 = hist([zeros(1,s_0_max) ones(1,s_1_max) 2*ones(1,s_2_max)],[0 1 2]);
h2_nominal_1 = h2_nominal_1/sum(h2_nominal_1);

h2_nominal_2 = hist([zeros(1,s_0_min) ones(1,s_1_min) 2*ones(1,s_2_min)],[0 1 2]);
h2_nominal_2 = h2_nominal_2/sum(h2_nominal_2);
% bar(h2, 'DisplayName', 'Time');
% legend('show');
% title('min-max')


bar([0 1 2],[h2_nominal_1;h2_nominal_2;h2_minmax]')
legend('show');