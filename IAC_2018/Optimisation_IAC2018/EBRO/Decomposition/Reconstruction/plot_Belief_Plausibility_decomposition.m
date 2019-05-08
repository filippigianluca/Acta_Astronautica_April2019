% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [LIST] = plot_Belief_Plausibility_decomposition(decomposition_end, in)


fmax    = [];
bpa_max = [];
fmin    = [];
bpa_min = [];

fmax_add = [];
bpa_add  = [];

position_FE_belief = [];
table = decomposition_end{1, 1}.final_FE;

for i = 1:length(decomposition_end)
    if isempty(decomposition_end{i})==0
        
        if in.flag_output.Belief     
            % BELIEF    
            
            for k = 1 : length(decomposition_end{i}.step_two)     % vector of all the f-values of the FEs                
                fmax    = [fmax     decomposition_end{i}.step_two{k}.F]; %decomposition_end{i}.step_two{k}.F;
                bpa_max = [bpa_max  decomposition_end{i}.step_two{k}.bpa]; 
                position_FE_belief = [position_FE_belief  decomposition_end{i}.step_two{k}.position]; 
                
            end    
                %%
%                 fmax_add = [fmax_add decomposition_end{i}.FE_add.F];
%                 bpa_add  = [bpa_add  decomposition_end{i}.FE_add.bpa];
        end
        
        if i ~= 1
            table = [table; decomposition_end{1, i}.final_FE(2:end,:)];
        end        
        
        if in.flag_output.Plausibility           
            % PLAUSIBILITY            
            for k = 1 : length(decomposition_end{i}.step_two)     % vector of all the f-values of the FEs                
                fmin    = [fmin     decomposition_end{i}.step_two{k}.F_Plausibility];
                bpa_min = [bpa_min  decomposition_end{i}.step_two{k}.bpa_Plausibility];
            end           
        end
                
    end
end



%  Belief
if in.flag_output.Belief 

% add the focal elements from an other threshold but with F<
    fmax    =  [fmax     fmax_add];
    bpa_max =  [bpa_max  bpa_add];
    
    %%
    [f_sorted_Bel, position_sorted_Bel] = sort(fmax);
    
    LIST.F_Bel    = [f_sorted_Bel(1)  f_sorted_Bel];
    LIST.Bel      = [0                cumsum(bpa_max(position_sorted_Bel))];
    LIST.position = position_FE_belief(position_sorted_Bel);
    LIST.table = table;
    %%

end


% Plausibility
if in.flag_output.Plausibility
    
    [f_sorted_Pl, position_sorted_Pl] = sort(fmin);
    
    LIST.F_Pl = [f_sorted_Pl(1)  f_sorted_Pl];
    LIST.Pl   = [0               cumsum(bpa_min(position_sorted_Pl))];
end





end