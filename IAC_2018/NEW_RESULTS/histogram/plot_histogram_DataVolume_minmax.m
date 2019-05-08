% history of data volume
% clear all; close all; clc

func_handle = @Acta_Astronautica_objconstr1;
threshold_vector = 600; %[500 600 700 800];
t = 0:365;

flag_plot_p2 = 0;
flag_plot_R2 = 0;
flag_plot_history =0;

ind = 1;


for threshold = threshold_vector
    clearvars -except threshold_vector t flag_plot_p2 flag_plot_R2 flag_plot_history ind threshold
    
    % load minmax_d and minmax_u
    load(strcat('Acta_nu',num2str(threshold),'_so_nfeval_minmax_3000000nfeval_inner_outer_20000'),'minmax');
    minmax_d = minmax.d;
    minmax_u = minmax.u{1};
    minmax_f = minmax.f{1};    
    
    par.fix.time = 365;
    TM = par.fix.time;
    
    
    
    [R1_function, R2_function, lam, mu] = new_histogram_Acta(minmax_d, minmax_u, par);
    

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
    

    
    
    %% probability distribution of R2
    
    s = mu+lam;
    p2 = mu/s + lam*exp(-t*s)/s;

    
    %%
    T = 0;
    T_fail_restore = 0;
%     V = R2;
    
    STATE = [];
    Time_state0 = [];
    Time_state1 = [];
    Time_state2 = [];
    for  kk = 1:500000


    
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
    h1 = hist(STATE,[0 1 2]);
    h1 = h1/sum(h1);
    bar(h1, 'DisplayName', 'N transitions');
    legend('show');
    title('min-max')
    
    figure
    sum_time_state0 = sum(Time_state0);
    sum_time_state1 = sum(Time_state1);
    sum_time_state2 = sum(Time_state2);
    
    s=[zeros(1,fix(sum_time_state0)) ones(1,fix(sum_time_state1)) 2.*ones(1,fix(sum_time_state2))];
    
    h2 = hist(s,[0 1 2]);
    h2 = h2/sum([sum_time_state0 sum_time_state1 sum_time_state2]);
    bar(h2, 'DisplayName', 'Time');
    legend('show');
    title('min-max')
    
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




