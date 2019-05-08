% history of data volume
% clear all; close all; clc



function [STATE, s] = histogram_nominal_2(rip_tot)


func_handle = @Acta_Astronautica_objconstr1;

t = 0:365;



ind = 1;


    clearvars -except threshold_vector t flag_plot_p2 flag_plot_R2 flag_plot_history ind threshold rip_tot
    
    % load minmax_d and minmax_u
    load('min_deterministic','d_min_det');
    load('worst_case_c_ddet_')

    d = d_min_det;
    u = worst_case.u{1, 1}  ;  
    
    par.fix.time = 365;
    TM = par.fix.time;
    
    
    
    [R1_function, R2_function, lam, mu] = new_histogram_Acta(d, u, par);
    

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

    
    %%
    T = 0;
    T_fail_restore = 0;
%     V = R2;
    
    STATE = [];
    Time_state0 = [];
    Time_state1 = [];
    Time_state2 = [];
    counter = 1;
    for  kk = 1:rip_tot
        
        a = [];
        b = []; 
        state = 2;
        
        state_2 = 1;
        state_1 = 0;
        T = 0;
        T_tot_fail = min(wblrnd(scale, beta));
        continua = 1;
        
        TM = 365;
        
        
        while T < min(TM,T_tot_fail) %&& continua 
            
%             T_partialfailure = 1e20;
%             T_restore = 0;

            
            %     while min(length(a), length(b)) <=1000
            
            aa = exprnd(1/lam); % lam: 2->1
            bb = exprnd(1/mu);  % mu:  1->2
            
            
            if bb<= aa
                a = [a aa];
                state = [2];
            else
                a = [a aa];
                b = [b bb-aa];
                state = [1 2];
            end


        
%                     continua = 0;
%                     if aa <= min((TM - T), (T_tot_fail-T)) && state_2
%                         a = [a aa];
%                         state_2 = 0;
%                         state_1 = 1;
%                         state = [state 1];
%                         continua = 1;
%                     end
%         
%                     if bb <= min((TM - T), (T_tot_fail-T)) && state_1
%                         b = [b bb];
%                         state_1 = 0;
%                         state_2 = 1;
%                         state = [state 2];
%                         continua = 1;
%                     end
            
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
            
            counter = counter + 1;
        end
        
        STATE = [STATE state];
        Time_state2 = [Time_state2 a];
        Time_state1 = [Time_state1 b];
        Time_state0 = [Time_state0 c];
    end
    
    
sum_time_state0 = sum(Time_state0);
sum_time_state1 = sum(Time_state1);
sum_time_state2 = sum(Time_state2);

s=[zeros(1,fix(sum_time_state0)) ones(1,fix(sum_time_state1)) 2.*ones(1,fix(sum_time_state2))];

return 