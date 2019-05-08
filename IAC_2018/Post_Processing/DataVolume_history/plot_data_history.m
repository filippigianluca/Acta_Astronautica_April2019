% history of data volume
% clear all; close all; clc


threshold_vector = 380;%[380 390 400 410 1000];
t = 0:365;

flag_plot_p2 = 0;
flag_plot_R2 = 0;
flag_plot_history =1;

ind = 1;


for threshold = threshold_vector
%% 
if threshold == 1000
    load('d_u_minmax_non_convergence.mat')
    d = d_minmax_non_convergence_1;
    u = u_minmax_non_convergence_1;
else
    load(strcat('minmax_nu',num2str(threshold)))
    d = dminmax;
    u = umaxC_dminmax;
    u = uminmax;
end
par.fix.time = 365;
TM = par.fix.time;



[R1_function, R2_function, lam, mu] = IAC2018_data_hystory(d, u, par);

R1= R1_function(100);
R2= R2_function(100);

if flag_plot_R2
    hold on
    p_R2(ind) = plot(t, R2_function(t));
end

%% weibul distribution for p0
beta = [ 0.7182
    0.3375
    1.4560
    0.356
    0.8874
    0.746
    0.5021
    0.4035
    0.3939];

scale = [ 3831
    6206945
    408
    21308746
    7983
    7733
    169.272
    1965868
    400982];
    
T_tot_fail = min(wblrnd(scale, beta));


%% probability distribution of R2

s = mu+lam;
p2 = mu/s + lam*exp(-t*s)/s;

if flag_plot_p2
    hold on
    p_p2(ind) = plot(t, p2,'DisplayName',strcat('\nu = ',num2str(threshold)));
end



%%
T = 0;
T_fail_restore = 0;
V = R2; 

while T < min(TM,T_tot_fail)
    
    T_partialfailure = 1e20;
    T_restore = 0;
    a = [];
    b = [];

%     while min(length(a), length(b)) <=1000 
    for  kk = 1:1000    
    aa = exprnd(1/lam);
    bb = exprnd(1/mu);
        if aa <=365
            a = [a aa];
        end
        if bb <= 365
            b = [b bb];
        end
    
    T_partialfailure = min(T_partialfailure, a(k));
    T_restore        = max(T_restore, b(k));
    end
     
    T_fail_restore = [T_fail_restore   sum(T_fail_restore)+T_partialfailure];
    T_fail_restore = [T_fail_restore   sum(T_fail_restore)+T_restore];
    
    V = [V R1 R2];
    T = sum(T_fail_restore);
end

if flag_plot_history
    hold on
    stairs(T_fail_restore, V)
end

%% plot total failure
% plot([T_tot_fail T_tot_fail], [0 R2], 'linewidth',2);

if flag_plot_p2
    lgd = legend(p_p2(1,:));
end
if flag_plot_R2
    lgd = legend(p_R2(1,:));
end

ind= ind+1;
end




