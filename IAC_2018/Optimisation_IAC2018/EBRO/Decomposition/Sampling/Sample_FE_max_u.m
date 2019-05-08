function [Sample] = Sample_FE_max_u(problem, minmax, minmin, coupled_position, num_coupled_vector, num_samples_tot, samples_input, Partial_curve, decomposition)


%% BELIEF
if problem.flag_output.Belief
    
    for n_sample = 1:length(samples_input)
        
        size_FEc = problem.num_interval(coupled_position);        
        if n_sample == 1
            pos_FEc_minmax  = minmax.FE.position_FE(coupled_position);
            
            
            bpa  = find_bpa(pos_FEc_minmax, problem.bpa{1}(coupled_position));
            n_FE = sub2ind_variable_dim(pos_FEc_minmax, size_FEc);     
            
            Sample{num_coupled_vector}.Belief_f            = 1;
            Sample{num_coupled_vector}.position            = [];
            Sample{num_coupled_vector}.f                   = minmax.f;
            Sample{num_coupled_vector}.num_coupled_vector  = num_coupled_vector;
            
            Sample{num_coupled_vector}.FE{1}.bpa       = bpa;    
            Sample{num_coupled_vector}.FE{1}.n_FE      = n_FE;     
            Sample{num_coupled_vector}.FE{1}.lb        = minmax.FE.lb(coupled_position);
            Sample{num_coupled_vector}.FE{1}.ub        = minmax.FE.ub(coupled_position);
            Sample{num_coupled_vector}.FE{1}.position  = minmax.FE.position_FE(coupled_position);
            Sample{num_coupled_vector}.FE{1}.upper_f   = minmax.f;
            Sample{num_coupled_vector}.FE{1}.upper_u   = minmax.u(coupled_position);
            
            
        else
            
          
            for i = 1:length(samples_input)
                ii = samples_input(i);
                FE_sample = problem.sample.FE_sample{ii}; 
                pos_FEc_sample  = problem.sample.FE_sample{ii}.position_FE(coupled_position);

            
                bpa  = find_bpa(pos_FEc_sample, problem.bpa{1}(coupled_position));
                n_FE = sub2ind_variable_dim(pos_FEc_sample, size_FEc);  
                
                pos_FE_PartialCurve = find(Partial_curve{1, 1}.Belief_FE_number  == n_FE);
                Sample{num_coupled_vector}.Belief_f(1+i)            = Partial_curve{num_coupled_vector}.Belief_FE_belief_partial(pos_FE_PartialCurve);
                Sample{num_coupled_vector}.position(1+i)           
                Sample{num_coupled_vector}.f(1+n)                   = decomposition{num_coupled_vector}.FocalElement{pos_FE_PartialCurve}.upper_f;

                Sample{num_coupled_vector}.FE{1}.bpa(1+i)       = bpa;    
                Sample{num_coupled_vector}.FE{1}.n_FE(1+i)      = n_FE;     
                Sample{num_coupled_vector}.FE{1}.lb{1+i}        = FE_sample.lb(coupled_position);
                Sample{num_coupled_vector}.FE{1}.ub{1+i}        = FE_sample.ub(coupled_position);
                Sample{num_coupled_vector}.FE{1}.position{1+i}  = FE_sample.position_FE(coupled_position);
                Sample{num_coupled_vector}.FE{1}.upper_f(1+i)   
                Sample{num_coupled_vector}.FE{1}.upper_u{1+i}               
            end
            
        end
        
        
        
    end
    
end




%% PLAUSIBILITY
if problem.flag_output.Plausibility
    
    %                 coupled_position = in.order_dim_u(sum(in.dim_u_i(1:in.num_functions)) + 1 : end);
    
    Sample{num_coupled_vector}.Plausibility_f     = 1;
    Sample{num_coupled_vector}.position           = [];
    Sample{num_coupled_vector}.f_Pl               = minmin.f;
    Sample{num_coupled_vector}.num_coupled_vector = num_coupled_vector;
    
    bpa = 1;
    pos_FE = '';
    for jj = coupled_position
        bpa = bpa*problem.bpa{1}{jj}(minmin.FE.position_FE(jj));
        pos_FE = strcat(pos_FE, ',',num2str(minmin.FE.position_FE(jj)));
    end
    Sample{num_coupled_vector}.FE_Plausibility{1}.bpa        = bpa;
    Sample{num_coupled_vector}.FE_Plausibility{1}.n_FE       = eval(strcat('sub2ind([',num2str(problem.num_interval(coupled_position)),']', pos_FE, ');'));
    Sample{num_coupled_vector}.FE_Plausibility{1}.lb         = minmin.FE.lb(coupled_position);
    Sample{num_coupled_vector}.FE_Plausibility{1}.ub         = minmin.FE.ub(coupled_position);
    Sample{num_coupled_vector}.FE_Plausibility{1}.position   = minmin.FE.position_FE(coupled_position);
    Sample{num_coupled_vector}.FE_Plausibility{1}.downer_f   = minmin.f;
    Sample{num_coupled_vector}.FE_Plausibility{1}.downer_u   = minmin.u(coupled_position);
end

return