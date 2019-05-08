function [Sample] = Sampling_Belief_Plausibility_coupled_u(Partial_curve, decomposition, problem, minmax, minmin, FE_minmax, FE_minmin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE each partial curve with the input-number of samples.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(strcat('Running Decomposition Sampling...'));


num_coupled_vector_tot = length(problem.dim_u_i(problem.num_functions+1:end));

for num_coupled_vector = 1:num_coupled_vector_tot
    
    num_vector_tot =  problem.num_functions + num_coupled_vector;              % number that identifies the coupled vector in all the vector (u1 u2 ... ui | u12 u13 ...uij)
    coupled_position = var2opt(num_vector_tot,problem);                        % components of vector u in the considered uij
    
    
    if problem.dim_u_i(num_vector_tot)>0
        
        num_samples_tot = problem.num_samples{num_coupled_vector};
        samples_input = problem.sample.FE_uSample_M{problem.indexing.matrix2linear{num_vector_tot}(1), problem.indexing.matrix2linear{num_vector_tot}(2)}; 
        num_samples_input = length(samples_input); 
        
%% new version        
%         if problem.flag_output.Belief
%             
%             % if the we take only one sample from the i-th partial curve, the
%             % sampling point will be the (min)-max solution
%             
% 
%             [Sample] = Sample_FE_max_u(problem, minmax, minmin, coupled_position, num_coupled_vector, num_samples_tot, samples_input, Partial_curve, decomposition);           
% 
% 
%             if  num_samples_input > 1
%                 [Sample] = Sample_FE_input_u(problem, minmax, minmin, coupled_position, num_coupled_vector);     
%             end
%             
%             if  num_samples_tot > num_samples_input 
%                 [Sample] = Sample_FE_input_u(problem, minmax, minmin, coupled_position, num_coupled_vector);     
%             end            
%         
%         end
%%         
        
  
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        if num_samples_tot == 1  %&& ~isempty(problem.input.u_sample{1})
            % BELIEF    
            if problem.flag_output.Belief
                
                
                
                Sample{num_coupled_vector}.Belief_f(1)           = 1;
                Sample{num_coupled_vector}.position(1)           = 1; % to change
                Sample{num_coupled_vector}.f(1)                  = minmax.f;               
                Sample{num_coupled_vector}.num_coupled_vector(1) = num_coupled_vector;
                
                bpa = 1;
                pos_FE = '';
                for jj = coupled_position
                    bpa = bpa*problem.bpa{1}{jj}(minmax.FE.position_FE(jj));
                    pos_FE = strcat(pos_FE, ',',num2str(minmax.FE.position_FE(jj)));
                end
                
                Sample{num_coupled_vector}.FE{1}.bpa(1)           = bpa; 
                
                if length(problem.num_interval(coupled_position)) == 1
                    Sample{num_coupled_vector}.FE{1}.n_FE  = minmax.FE.position_FE(jj);
                else
                    Sample{num_coupled_vector}.FE{1}.n_FE  = eval(strcat('sub2ind([',num2str(problem.num_interval(coupled_position)),']', pos_FE, ');'));
                end
                
                Sample{num_coupled_vector}.FE{1}.lb        = minmax.FE.lb(coupled_position);
                Sample{num_coupled_vector}.FE{1}.ub        = minmax.FE.ub(coupled_position);
                Sample{num_coupled_vector}.FE{1}.position  = minmax.FE.position_FE(coupled_position);
                Sample{num_coupled_vector}.FE{1}.upper_f   = minmax.f;
                Sample{num_coupled_vector}.FE{1}.upper_u   = minmax.u(coupled_position);
  
                
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                % if sampling points in input in a different FE than the
                % minmax
%                 ind_sample = 1;
%                 if ~all(problem.FE_minmax_uc(:,num_coupled_vector))
%                     for kk = 1:length(problem.N_ind_components_u{num_vector_tot})       %        size(problem.FE_minmax_uc,1)
%                         
%                         kkk = problem.N_ind_components_u{num_vector_tot}(kk);
%                         if problem.FE_minmax_uc(kkk,num_coupled_vector) == 0
%                             
%                             
%                             pos_FE = '';
%                             for jj = coupled_position
%                                 pos_FE = strcat(pos_FE, ',',num2str(problem.FE_sample{kkk}.position_FE(jj)));
%                             end
%                             
%                             
%                             Sample{num_coupled_vector}.FE{1}.bpa              = bpa; 
%                             if length(problem.num_interval(coupled_position)) == 1
%                                 Sample{num_coupled_vector}.FE{1+ind_sample}.n_FE  = pos_FE;
%                             else
%                                 Sample{num_coupled_vector}.FE{1+ind_sample}.n_FE  = eval(strcat('sub2ind([',num2str(problem.num_interval(coupled_position)),']', pos_FE, ');'));
%                             end                            
%                             pos_FE_int = Sample{num_coupled_vector}.FE{1+ind_sample}.n_FE;
%                             
%                             
%                             Sample{num_coupled_vector}.Belief_f(1+ind_sample)           = Partial_curve{num_coupled_vector}.Belief_FE_belief_partial(pos_FE_int);
% %                             Sample{num_coupled_vector}.position(1+ind_sample)           = decomposition{num_coupled_vector}.FocalElement{pos_FE_int}.position; 
%                             Sample{num_coupled_vector}.f(1+ind_sample)                  = decomposition{num_coupled_vector}.FocalElement{pos_FE_int}.upper_f;               
% 
%                             
%                             Sample{num_coupled_vector}.FE{1+ind_sample}  =  decomposition{num_coupled_vector}.FocalElement{pos_FE_int};
%                           
%                             
%                             
%                         end
%                         ind_sample = ind_sample + 1;
%                     end
%                     
%                 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                

            end
            
            % PLAUSIBILITY    
            if problem.flag_output.Plausibility
                
%                 coupled_position = in.order_dim_u(sum(in.dim_u_i(1:in.num_functions)) + 1 : end);
                
                Sample{num_coupled_vector}.Plausibility_f =  1;
                Sample{num_coupled_vector}.position = 1; % to change
                Sample{num_coupled_vector}.f_Pl = minmin.f;               
                Sample{num_coupled_vector}.num_coupled_vector = num_coupled_vector;
                
                bpa = 1;
                pos_FE = '';
                for jj = coupled_position
                    bpa = bpa*problem.bpa{1}{jj}(minmin.FE.position_FE(jj));
                    pos_FE = strcat(pos_FE, ',',num2str(minmin.FE.position_FE(jj)));
                end
                Sample{num_coupled_vector}.FE_Plausibility{1}.bpa       = bpa; 
                Sample{num_coupled_vector}.FE_Plausibility{1}.n_FE      = eval(strcat('sub2ind([',num2str(problem.num_interval(coupled_position)),']', pos_FE, ');'));
                Sample{num_coupled_vector}.FE_Plausibility{1}.lb        = minmin.FE.lb(coupled_position);
                Sample{num_coupled_vector}.FE_Plausibility{1}.ub        = minmin.FE.ub(coupled_position);
                Sample{num_coupled_vector}.FE_Plausibility{1}.position  = minmin.FE.position_FE(coupled_position);
                Sample{num_coupled_vector}.FE_Plausibility{1}.downer_f   = minmin.f;
                Sample{num_coupled_vector}.FE_Plausibility{1}.downer_u   = minmin.u(coupled_position);
            end
            
            
            
            
             
            
            
        else            
            %% BELIEF           
            if problem.flag_output.Belief
                
                
                [Bel_unique_values, index] = sort(Partial_curve{num_coupled_vector}.Belief_FE_belief_partial(:));
                
                
                if num_samples_tot > length(Bel_unique_values)
                    
                    num_samples_tot = length(Bel_unique_values);
                end
                
                Bel_subset_num = floor((length(Bel_unique_values)-1)/(num_samples_tot-1));
                
                
                
                
                
                for i = 1:num_samples_tot
                    flag_continue = 1;
                    N_Focal_Element_coupled_vector_k = length(decomposition{1,num_coupled_vector}.FocalElement);
                    
                    if i==1
                        Sample{num_coupled_vector}.Belief_f(num_samples_tot+1-i) =  1;
                        fun(i) = Partial_curve{num_coupled_vector}.Belief_FE_function(index(end));
                        Sample{num_coupled_vector}.position(num_samples_tot+1-i) = N_Focal_Element_coupled_vector_k;
                    else
                        
                        
                        
                        %                     Bel_subset = Bel_unique_values(end-(i-1)*Bel_subset_num : end-(i-2)*Bel_subset_num -1);
                        %                     [Sample{num_coupled_vector}.Belief_f(num_samples_tot+1-i), idx] = datasample(Bel_subset, 1);
                        %
                        %                     fun(i) = Partial_curve{num_coupled_vector}.Belief_FE_function(index(end-(i-1)*Bel_subset_num + idx -1));
                        %
                        %                     Sample{num_coupled_vector}.position(num_samples_tot+1-i) = index(end-(i-1)*Bel_subset_num + idx -1);% - par_treshold; par_treshold;% index(par_treshold);%
                        % global par_treshold
                        
                        % the following two lines make the samples from the right to the left and with step =1 over the focal elements.
%                         flag_continue = 1;
                        if isempty(problem.sample.FE_sample{1, 1})
                            Sample{num_coupled_vector}.Belief_f(num_samples_tot+1-i) = Partial_curve{num_coupled_vector}.Belief_FE_belief_partial(index(end-i+1)); %(index(par_treshold));%
                            fun(i) = Partial_curve{num_coupled_vector}.Belief_FE_function(index(end-i+1)); %(index(par_treshold));% sample in the right
                            Sample{num_coupled_vector}.position(num_samples_tot+1-i) = index(end-i+1);%index(par_treshold);%
                        else
                            pos_c_u_sample_matrix = problem.indexing.matrix2linear{num_vector_tot};
                            pos_c_u_sample = problem.CIM{pos_c_u_sample_matrix(1), pos_c_u_sample_matrix(2)};
                            FE_pos_c_u_sample = problem.sample.FE_sample{1, 1}.position_FE(pos_c_u_sample);
                            ind = sub2ind_variable_dim(FE_pos_c_u_sample, 2*ones(1,length(FE_pos_c_u_sample)));
                            
%                             if ~any(Sample{num_coupled_vector}.position == ind)
                                Sample{num_coupled_vector}.Belief_f(num_samples_tot+1-i) = Partial_curve{num_coupled_vector}.Belief_FE_belief_partial(ind);
                                fun(i)                                                   = Partial_curve{num_coupled_vector}.Belief_FE_function(ind);           
                                Sample{num_coupled_vector}.position(num_samples_tot+1-i) = index(ind);   
%                             else
%                                 flag_continue = 0;
                            
%                             end
                        end
                    end
                    
                    
                    
                    
%                     if flag_continue
                    
                        for j = 1 : N_Focal_Element_coupled_vector_k

                            % Belief
                            if abs(decomposition{1,num_coupled_vector}.FocalElement{j}.upper_f - fun(i)) == 0

                                Sample{num_coupled_vector}.FE{num_samples_tot+1-i} = decomposition{1,num_coupled_vector}.FocalElement{j};

                            end


                        end
                    
%                     end
                end
                
%                 Sample{num_coupled_vector}.FE       = Sample{num_coupled_vector}.FE(~cellfun('isempty',Sample{num_coupled_vector}.FE));            
%                 Sample{num_coupled_vector}.Belief_f = Sample{num_coupled_vector}.Belief_f(Sample{num_coupled_vector}.Belief_f~=0);    
%                 Sample{num_coupled_vector}.position = Sample{num_coupled_vector}.position(Sample{num_coupled_vector}.position~=0);
                
                Sample{num_coupled_vector}.f = fliplr(fun);      
                Sample{num_coupled_vector}.num_coupled_vector = num_coupled_vector;
              
                
            end
            % par_treshold = par_treshold +1;
            
            
            
            %% PLAUSIBILITY
            if problem.flag_output.Plausibility
                
                
                
                [Plausibility_unique_values, index_Pl] = sort(Partial_curve{num_coupled_vector}.Plausibility_FE_plausibility_partial(:));
                
                if num_samples_tot >  length(Plausibility_unique_values)
                    
                    num_samples_tot = length(Plausibility_unique_values);
                end
                
                Plausibility_subset_num = floor((length(Plausibility_unique_values)-1)/(num_samples_tot-1));
                
                
                
                for i = 1:num_samples_tot
                    
                    if i==1
                        
                        Sample{num_coupled_vector}.Plausibility_f(num_samples_tot+1-i) =  1;
                        fun_Pl(i) = Partial_curve{num_coupled_vector}.Plausibility_FE_function(index_Pl(end));
                    else
                        
                        Plausibility_subset = Plausibility_unique_values(end-(i-1)*Plausibility_subset_num : end-(i-2)*Plausibility_subset_num -1);
                        %                     [Sample{num_coupled_vector}.Plausibility_f(num_samples_tot+1-i), idx_Pl] = datasample(Plausibility_subset, 1);
                        %
                        %                     fun_Pl(i) = Partial_curve{num_coupled_vector}.Plausibility_FE_function(index_Pl(end-(i-1)*Plausibility_subset_num + idx_Pl -1));
                        
                        % the following two lines make the samples from the right to the left and with step =1 over the focal elements.
                        Sample{num_coupled_vector}.Plausibility_f(num_samples_tot+1-i)=Partial_curve{num_coupled_vector}.Plausibility_FE_plausibility_partial(index_Pl(end-i+1));
                        fun_Pl(i) = Partial_curve{num_coupled_vector}.Plausibility_FE_function(index_Pl(end-i+1)); % sample in the right
                    end
                    
                    for j=1:length(decomposition{1,num_coupled_vector}.FocalElement)
                        
                        if abs(decomposition{1,num_coupled_vector}.FocalElement{j}.downer_f - fun_Pl(i)) == 0 %< (minmax.f-minmin.f)/1000
                            
                            Sample{num_coupled_vector}.FE_Plausibility{num_samples_tot+1-i} = decomposition{1,num_coupled_vector}.FocalElement{j};
                            
                        end
                        
                    end
                    
                    
                end
                
                Sample{num_coupled_vector}.f_Pl = fun_Pl;
                
                Sample{num_coupled_vector}.num_coupled_vector = num_coupled_vector;
                
            end
            
        end
        
        
        
        
        
    end
    
end

end