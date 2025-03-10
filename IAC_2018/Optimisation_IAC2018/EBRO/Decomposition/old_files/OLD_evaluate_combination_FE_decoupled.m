function [decomposition_end] = OLD_evaluate_combination_FE_decoupled(decomposition_end, num_sample, problem_decomposition, minmax, minmin, u_max_tot, Sample, PS, u_max_tot_Plausibility, n_obj, n_point, in, Partial_curve)


%% FIND number of combinations for the final curve
num_combinations_max = 1;
num_functions_tot = length(decomposition_end{num_sample}.step_one);

num_functions_tot_new = 0;
for num_functions = 1:num_functions_tot
    if isempty(decomposition_end{num_sample}.step_one{num_functions})==0
        num_combinations_max = num_combinations_max*length(decomposition_end{num_sample}.step_one{num_functions});
        
        num_functions_tot_new = num_functions_tot_new + 1;
        sub_not_empty(num_functions_tot_new) = num_functions;
    end
end
%num_functions_tot = num_functions_tot_new;

%% BELIEF
decomposition_end{num_sample}.final_FE = [];

if in.output == 0 || in.output == 2
    
    d_belief = minmax.d;
    f_cost = problem_decomposition.objfun{n_obj}(d_belief, u_max_tot, problem_decomposition);
    
    
    
    % EVALUATE bpa-scale: for each sample a scaled curve is evaluated.
    
    bpa_scale = 1;
    Sample_u_lb = [];
    Sample_u_ub = [];
    
    for j = 1:length(PS)
        if in.dim_u(in.num_functions + j)>0
            if PS(j) == 1
                bpa_scale = bpa_scale*Sample{j}.Belief_f(PS(j));
                
            else
                bpa_scale = bpa_scale*(Sample{j}.Belief_f(PS(j)) - Sample{j}.Belief_f(PS(j)-1));
                
            end
            Sample_u_lb = [Sample_u_lb  Sample{j}.FE{1, 1}.lb];
            Sample_u_ub = [Sample_u_ub  Sample{j}.FE{1, 1}.ub];
        end
    end
    
    
    % RECONSTRUCTION of all the focal elements of the final scaled curves
   
    % initialise for add-FE
%     bpa_add = 0;
    F_max = 0;
    F_ij_add2Delta_ij = [];
    bpa_add2Delta_ij  = [];
    N_add = 1;
    
    for num_combination = 1:num_combinations_max
        
        [P_FE] = position_sum_function(num_combinations_max, num_functions_tot, num_combination, decomposition_end, num_sample);
        
        bpa =1;
        F = 0;
        upper_u_tot_decoupled = [];
        lower_bound_tot_decoupled = [];
        upper_bound_tot_decoupled = [];
        
        
        for i = 1:length(P_FE)
            bpa = bpa*decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.bpa;
            F = F + decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.upper_f;
            upper_u_tot_decoupled = [upper_u_tot_decoupled decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.upper_u];
            lower_bound_tot_decoupled = [lower_bound_tot_decoupled decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.lb];
            upper_bound_tot_decoupled = [upper_bound_tot_decoupled decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.ub];
        end
        lower_bound = [lower_bound_tot_decoupled  Sample_u_lb]; %Sample{num_sample}.FE{1, 1}.lb];
        upper_bound = [upper_bound_tot_decoupled  Sample_u_ub]; %Sample{num_sample}.FE{1, 1}.ub];        
        
        N=0;
        for i=1:num_functions_tot
            if in.dim_u(i)>0
                N = N+1;
            end
        end
        decomposition_end{num_sample}.step_two{num_combination}.F = F - (N-1)*f_cost;
        
        F_belief(num_combination) = decomposition_end{num_sample}.step_two{num_combination}.F;
        
        if F_belief(num_combination) > F_max
            position_F_max = num_combination;
            F_max = decomposition_end{num_sample}.step_two{num_combination}.F;
        end
        
        P_partial = [];
        %         for ii=1:length(decomposition_end{num_sample}.step_one)
        %               P_partial = [P_partial decomposition_end{num_sample}.step_one{ii}{1, 1}.n_FE  Sample{num_sample}.FE{num_sample}];
    for j = 1:length(PS)
        if in.dim_u(in.num_functions + j)>0
        P_partial = [P_partial Partial_curve{j}.Belief_FE_number(PS(j))];
        end
    end
        %         end
        
   
        %% change (?) the 2 "if" condition
        coupled_vectors = in.dim_u(in.num_functions +1 : end);
        position_non_zero_c_v = find(coupled_vectors>0);
        num_coupled_vectors = length(coupled_vectors);%(coupled_vectors>0));
        num_comb_add = 1;
        
        for jj = position_non_zero_c_v
            num_comb_add = num_comb_add*length(Partial_curve{jj}.Belief_FE_function);
        end
  %% 
    position = [];
    for ii=1:length(Sample)
        if isempty(Sample{ii}) == 0
            position = [position Sample{ii}.position(PS(ii))]; 
        end
    end
    
    position_next = position;
    for ii=1:length(Sample)
        if isempty(Sample{ii}) == 0
            if Sample{ii}.position(PS(ii))~=Sample{ii}.position(end)
                position_next = [position_next Sample{ii}.position(PS(ii)+1)]; 
            end
        end
    end
    
    mm_position = 1;
    for mm=1:length(Partial_curve)
        if isempty(Partial_curve{mm}) == 0
            structure{mm}.FE = Partial_curve{mm}.Belief_FE_function;
            mm_position = mm_position +1;
        end
    end
    for num_sample_add = 1 : num_comb_add
        [PS_add] = position_Sample(num_comb_add, num_coupled_vectors, num_sample_add, structure, in);
 
        
        if sum(PS_add(PS_add>0) > position)==1% || all(PS_add(PS_add>0) <= position)==0% && min(PS_add(PS_add>0) < position_next)==1   %Sample{num_sample}.position(position_non_zero_c_v))==1
            coupled_add = [];
            Delta_belief_add = 1;
            PS_add
            for kk = 1: length(PS_add)
                 if isempty(Partial_curve{kk}) == 0
                    coupled_add = [coupled_add  Partial_curve{kk}.u_coupled(PS_add(kk),:)];
                    Delta_belief_add = Delta_belief_add*Partial_curve{kk}.Belief_FE_belief(PS_add(kk));
                 end
            end
            
            u_ij = [upper_u_tot_decoupled coupled_add];
            F_ij = problem_decomposition.objfun{n_obj}(d_belief, u_ij, problem_decomposition);
            
         if F_ij < Sample{num_sample}.f

                    F_ij_add2Delta_ij(N_add) = F_ij;
                    bpa_add2Delta_ij(N_add)  = bpa*Delta_belief_add;

%                   bpa_add = bpa_add + bpa*Partial_curve{1, 1}.Belief_FE_belief(N_ij);
                    N_add = N_add + 1;
         end
        end
    end

    
    
coupled = [];
for kkk = 1: length(PS)
    if in.dim_u(in.num_functions + kkk)>0
        coupled = [coupled  Partial_curve{kkk}.u_coupled(PS(kkk),:)];
%         Delta_belief_add = Delta_belief_add*Partial_curve{kkk}.Belief_FE_belief(PS_add(kkk));
    end
end    
 
% OLD
% %         if decomposition_end{num_sample}.step_two{num_combination}.F < Sample{num_sample}.f
%             for N_ij = Sample{j}.position + 1 : length(Partial_curve{j}.Belief_FE_belief_partial)
%                 
%                 u_ij = [upper_u_tot_decoupled Partial_curve{j}.u_coupled(N_ij,:)];
%                 
%                 F_ij = problem_decomposition.objfun{n_obj}(d_belief, u_ij, problem_decomposition);
% %                 if F_ij < Sample{num_sample}.f
%                     F_ij_add2Delta_ij(N_add) = F_ij;
%                     bpa_add2Delta_ij(N_add)  = bpa*Partial_curve{j}.Belief_FE_belief(N_ij);
% %                     bpa_add = bpa_add + bpa*Partial_curve{1, 1}.Belief_FE_belief(N_ij);
%                     N_add = N_add + 1;
% %                 end
%             end
% %         end
bpa_2_remuve = 0;
% if num_sample>1
%      bpa_2_remuve = decomposition_end{num_sample-1}.FE_add.bpa;
% end
        decomposition_end{num_sample}.step_two{num_combination}.bpa = bpa*bpa_scale - bpa_2_remuve;
        
        final_FE = [P_FE P_partial];
        decomposition_end{num_sample}.step_two{num_combination}.FE = final_FE;
        decomposition_end{num_sample}.step_two{num_combination}.F = F - (num_functions_tot-1)*f_cost;
        decomposition_end{num_sample}.step_two{num_combination}.u_max = [upper_u_tot_decoupled coupled];%Partial_curve{1}.u_coupled(Sample{1, 1}.position(num_sample),:)];
  
        
        decomposition_end{num_sample}.final_FE = [decomposition_end{num_sample}.final_FE;...
                                                  {final_FE}...
                                                  F_belief(num_combination)...
                                                  {decomposition_end{num_sample}.step_two{num_combination}.u_max}...
                                                  decomposition_end{num_sample}.step_two{num_combination}.bpa...
                                                  {lower_bound}...
                                                  {upper_bound}];        
        
    end
 
    
    % add bpa contribute of the FE not considered with the threshold
    if isempty(F_ij_add2Delta_ij)==0
%         [F_add_sort, position_F_add_sort] = sort(F_ij_add2Delta_ij);
%         [F_sort,     position_F_sort]     = sort(F_belief);        
%         for k = 1 : num_combinations_max
%             
%             if k ==1
%                 bpa_to_sum = find(F_ij_add2Delta_ij<F_sort(k));
%             else
%                 bpa_to_sum = find(F_ij_add2Delta_ij<F_sort(k) & F_ij_add2Delta_ij>=F_sort(k-1)==1);
%                 
%             end
%             
%             bpa_add_end(position_F_sort(k)) = sum(bpa_add2Delta_ij(bpa_to_sum));
%             
%             %         decomposition_end{num_sample}.step_two{position_F_max}.bpa = decomposition_end{num_sample}.step_two{position_F_max}.bpa + bpa_add;
%             
%         end
        
%         for ll = 1 : num_combinations_max
%             decomposition_end{num_sample}.step_two{num_combination}.bpa = decomposition_end{num_sample}.step_two{num_combination}.bpa + bpa_add_end(ll);
%         end
    end
    
for i=1:length(decomposition_end{1, 1}.step_two)
    a(i)=decomposition_end{1, 1}.step_two{i}.F;
end


[aa, pos_aa] = find(F_ij_add2Delta_ij>=min(a));
[aaa, pos_aaa] = find(F_ij_add2Delta_ij==min(a));
[bb, pos_bb] = find(F_ij_add2Delta_ij<min(a));

bpa_bb = sum(bpa_add2Delta_ij(pos_bb));
bpa_add2Delta_ij(pos_aaa) = bpa_add2Delta_ij(pos_aaa(1)) + bpa_bb;

decomposition_end{num_sample}.FE_add.F = aa;                     % F_ij_add2Delta_ij
decomposition_end{num_sample}.FE_add.bpa = bpa_add2Delta_ij(aa); % bpa_add2Delta_ij



end

%% PLAUSIBILITY

if in.output == 1 || in.output == 2
    
    d_plausibility = minmin.d;
    f_cost_Plausibility = problem_decomposition.objfun{n_obj}(d_plausibility, u_max_tot_Plausibility, problem_decomposition);
    
    
    % EVALUATE bpa-scale: for each sample a scaled curve is evaluated.
    
    bpa_scale_Plausibility = 1;
    for j = 1:length(PS)
        if in.dim_u(in.num_functions + j)>0
            if PS(j) == 1
                
                bpa_scale_Plausibility = bpa_scale_Plausibility*Sample{j}.Plausibility_f(PS(j));
            else
                
                bpa_scale_Plausibility = bpa_scale_Plausibility*(Sample{j}.Plausibility_f(PS(j)) - Sample{j}.Plausibility_f(PS(j)-1));
            end
        end
    end
    
    
    % RECONSTRUCTION of all the focal elements of the final scaled curves
    
    for num_combination = 1:num_combinations_max
        
        [P_FE] = position_sum_function(num_combinations_max, num_functions_tot, num_combination, decomposition_end, num_sample);
        
        bpa =1;
        
        F_Plausibility = 0;
        for i = 1:length(P_FE)
            bpa = bpa*decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.bpa;
            
            F_Plausibility = F_Plausibility + decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.downer_f;
        end
        
        
        
        decomposition_end{num_sample}.step_two{num_combination}.bpa_Plausibility = bpa*bpa_scale_Plausibility;
        
        
        
        decomposition_end{num_sample}.step_two{num_combination}.F_Plausibility = F_Plausibility - (num_functions_tot-1)*f_cost_Plausibility;
    end
    
    
end



end