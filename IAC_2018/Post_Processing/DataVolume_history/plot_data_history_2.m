% history of data volume
% clear all; close all; clc


threshold_vector = [380 390 400 410];
t = 0:365;

flag_plot_p2 = 0;
flag_plot_R2 = 0;
flag_plot_history =0;

ind = 1;


for threshold = threshold_vector
    clearvars -except threshold_vector t flag_plot_p2 flag_plot_R2 flag_plot_history ind threshold
    
    %%
    if threshold == 1000
        load('d_u_minmax_non_convergence.mat')
        d = d_minmax_non_convergence_1;
        u = u_minmax_non_convergence_1;
    elseif threshold == 2000
        load('minimum_d_nominal_u')
        load('ebro_IAC_max_CONSTR_dmin_nominal_u','worst_case')
        d = d_min_u_nominal;
        u = worst_case.u{1, 1};        
    else
        load(strcat('minmax_nu',num2str(threshold)))
        d = dminmax;
%         u = umaxC_dminmax;
        u = uminmax;
    end
    par.fix.time = 365;
    TM = par.fix.time;
    
    
    
    [R1_function, R2_function, lam, mu] = IAC2018_data_hystory(d, u, par);
    
%     R1= R1_function(100);
%     R2= R2_function(100);
    
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
%     V = R2;
    
    STATE = [];
    Time_state0 = [];
    Time_state1 = [];
    Time_state2 = [];
    for  kk = 1:100000
        
        a = [];
        b = []; 
        state = 2;
        
        state_2 = 1;
        state_1 = 0;
        T = 0;
        T_tot_fail = min(wblrnd(scale, beta));
        continua = 1;
        
        TM = 365;
        
        while T < min(TM,T_tot_fail) && continua 
            
%             T_partialfailure = 1e20;
%             T_restore = 0;

            
            %     while min(length(a), length(b)) <=1000
            
            aa = exprnd(1/lam); % lam: 2->1
            bb = exprnd(1/mu);  % mu:  1->2
            
            continua = 0;
            if aa <= min((TM - T), (T_tot_fail-T)) && state_2
                a = [a aa];
                state_2 = 0;
                state_1 = 1;
                state = [state 1];
                continua = 1;
            end
            
            if bb <= min((TM - T), (T_tot_fail-T)) && state_1
                b = [b bb];
                state_1 = 0;
                state_2 = 1;
                state = [state 2];
                continua = 1;
            end
            
%             T_partialfailure = min(T_partialfailure, a(end);
%             T_restore        = max(T_restore, b);
%             
%             
%             T_fail_restore = [T_fail_restore   sum(T_fail_restore)+T_partialfailure];
%             T_fail_restore = [T_fail_restore   sum(T_fail_restore)+T_restore];
%             
%             V = [V R1 R2];
            T = sum(a)+sum(b);
            c=[];
            if T>=T_tot_fail && T<=TM
                state = [state 0];
                c = TM-T_tot_fail;
                break
            end
            
            if T<T_tot_fail && T<=TM && ~continua
                
                if state_1
                    b = [b T_tot_fail-T];
                end
                if state_2
                    a = [a T_tot_fail-T];
                end
                break
            end
        end
        
        STATE = [STATE state];
        Time_state2 = [Time_state2 a];
        Time_state1 = [Time_state1 b];
        Time_state0 = [Time_state0 c];
    end
    
    figure
    hist(STATE,[0 1 2]);
    title(['N nu',num2str(threshold)])
    
    figure
    sum_time_state0 = sum(Time_state0);
    sum_time_state1 = sum(Time_state1);
    sum_time_state2 = sum(Time_state2);
    
    s=[zeros(1,fix(sum_time_state0)) ones(1,fix(sum_time_state1)) 2.*ones(1,fix(sum_time_state2))];
    hist(s,[0 1 2]);
    title(['time nu',num2str(threshold)])
    
    
    
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




