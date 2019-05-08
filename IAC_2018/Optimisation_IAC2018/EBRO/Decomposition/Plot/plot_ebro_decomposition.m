function [] = plot_ebro_decomposition(problem, Partial_curve, LIST, LIST_EXACT)

hold on


ind_plot = 1;

for i = problem.num_functions +1 : length(problem.dim_u_i)
    if problem.dim_u_i(i)>0 && problem.num_samples{i-problem.num_functions}  ~= 1
        
        % partial Belief curves with Decomposition approach
        if problem.flag_output.Belief           
            Bel_start = Partial_curve{i - problem.num_functions}.Belief_FE_function(1);
            stairs([Bel_start; Partial_curve{i - problem.num_functions}.Belief_FE_function], [0 Partial_curve{i - problem.num_functions}.Belief_FE_belief_partial],'b','linewidth',1);
            legend('partial Belief (ENM)')
        end
        
        % partial Plausibility curves with Decomposition approach
        if problem.flag_output.Plausibility           
            Pl_start = Partial_curve{i - problem.num_functions}.Plausibility_FE_function(1);
            stairs([Pl_start; Partial_curve{i - problem.num_functions}.Plausibility_FE_function], [0 Partial_curve{i - problem.num_functions}.Plausibility_FE_plausibility_partial],'r', 'linewidth',1)  
        end
    end
end





% final Belief curve with Decomposition approach
if problem.flag_output.Belief
    s(ind_plot) = stairs(LIST.F_Bel, LIST.Bel,'b', 'linewidth',2,'DisplayName','Belief (ENM)');
    ind_plot = ind_plot + 1;
end


% final Plausibility curve with Decomposition approach
if problem.flag_output.Plausibility
    s(ind_plot) = stairs(LIST.F_Pl, LIST.Pl,'r','linewidth',2,'DisplayName','Plausibility (ENM)');
    ind_plot = ind_plot + 1;
end


% exact Belief curve
if problem.flag_output.exact_Belief
    s(ind_plot) = stairs(LIST_EXACT.F_Bel, LIST_EXACT.Bel, 'k', 'linewidth', 1,'DisplayName','Belief (exact)');
    ind_plot = ind_plot + 1;
end


% exact Plausibility curve
if problem.flag_output.exact_Plausibility
    s(ind_plot) = stairs(LIST_EXACT.F_Pl, LIST_EXACT.Pl, 'k', 'linewidth',1,'DisplayName','Plausibility (exact)');
end



legend(s(:));

hold off













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot of the partial curves
% for i = problem.num_functions +1 : length(problem.dim_u_i)
%     
%     hold on
%     if problem.dim_u_i(i)>0
%         
%         if problem.flag_output.Belief
%             Bel_start = Partial_curve{i - problem.num_functions}.Belief_FE_function(1);
%             stairs([Bel_start; Partial_curve{i - problem.num_functions}.Belief_FE_function], [0 Partial_curve{i - problem.num_functions}.Belief_FE_belief_partial],'b','linewidth',1)
%             xlabel('F')
%             ylabel('partial Belief')
%             title('Exchange variables (decomposition approach)')
%         elseif problem.flag_output.Plausibility
%             Pl_start = Partial_curve{i - problem.num_functions}.Plausibility_FE_function(1);
%             stairs([Pl_start; Partial_curve{i - problem.num_functions}.Plausibility_FE_function], [0 Partial_curve{i - problem.num_functions}.Plausibility_FE_plausibility_partial],'r', 'linewidth',1)
%             xlabel('F')
%             ylabel('partial Plausibility')
%             title('Exchange variables (decomposition approach)')
%         elseif problem.flag_output.Belief && problem.flag_output.Plausibility
%             Bel_start = Partial_curve{i - problem.num_functions}.Belief_FE_function(1);
%             Pl_start = Partial_curve{i - problem.num_functions}.Plausibility_FE_function(1);
%             stairs([Bel_start; Partial_curve{i - problem.num_functions}.Belief_FE_function], [0 Partial_curve{i - problem.num_functions}.Belief_FE_belief_partial],'b', 'linewidth',1)
%             hold on
%             stairs([Pl_start; Partial_curve{i - problem.num_functions}.Plausibility_FE_function], [0 Partial_curve{i - problem.num_functions}.Plausibility_FE_plausibility_partial],'r', 'linewidth',1)
%             xlabel('F')
%             ylabel('partial Belief & Plausibility')
%             title('Exchange variables (decomposition approach)')
%         end
%         
%     end
%     
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             figure
%             hold on
%             if problem.output == 0
%                 stairs(LIST_EXACT.F_Bel, LIST_EXACT.Bel, 'k', 'linewidth', 1)
%                 stairs(LIST.F_Bel, LIST.Bel,'b', 'linewidth',2)
%                 legend('Belief Exact','Belief Decomposition')
%             elseif problem.output == 1
%                 stairs(LIST_EXACT.F_Pl, LIST_EXACT.Pl, 'k', 'linewidth',1)
%                 stairs(LIST.F_Pl, LIST.Pl,'r','linewidth',2)
%                 legend('Plausibility Exact','Plausibility Decomposition')
%             elseif problem.output == 2
%                 stairs(LIST_EXACT.F_Bel, LIST_EXACT.Bel, 'k', 'linewidth', 1)
%                 stairs(LIST_EXACT.F_Pl, LIST_EXACT.Pl, 'k', 'linewidth', 1)
%                 stairs(LIST.F_Bel, LIST.Bel,'b','linewidth',1)
%                 stairs(LIST.F_Pl, LIST.Pl,'r','linewidth',1)
%                 legend('Belief e Plausibility Exact','Belief e Plausibility Decomposition')
%             end
%             hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             hold on
%             if problem.output == 0
%                 stairs(LIST_EXACT.F_Bel, LIST_EXACT.Bel, 'b', 'linewidth', 3)
%                 legend('Belief Exact','Belief Decomposition')
%             elseif problem.output == 1
%                 stairs(LIST_EXACT.F_Pl, LIST_EXACT.Pl, 'r', 'linewidth', 3)
%                 legend('Plausibility Exact','Plausibility Decomposition')
%             elseif problem.output == 2
%                 stairs(LIST_EXACT.F_Bel, LIST_EXACT.Bel, 'b', 'linewidth', 3)
%                 stairs(LIST_EXACT.F_Pl, LIST_EXACT.Pl, 'r', 'linewidth', 3)
%                 legend('Belief e Plausibility Exact','Belief e Plausibility Decomposition')
%             end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




return